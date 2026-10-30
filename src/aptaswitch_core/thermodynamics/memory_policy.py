"""Admission control for native folding workspaces, independent of the model.

Available-memory observations are advisory, never a promise against allocations
by other processes. Reserve one quarter (at least 256 MiB) for the host, and
account for concurrent requests in this process before allocating native data.
The native engine independently enforces the per-request byte limit.
"""
from __future__ import annotations

from contextlib import contextmanager
import ctypes as C
from functools import lru_cache
import os
from pathlib import Path, PurePosixPath
import sys
import threading

_MIB = 1024**2
_UNKNOWN_BUDGET = 512 * _MIB
_LOCK = threading.Lock()
_RESERVED = 0


def workspace_bytes(n: int, *, tables: int = 3) -> int:
    """Conservative connected-DP estimate for the requested observable.

    Each table has twenty-one triangular fields. Partition and sampling need
    one table; pairs need two. Eager builds' minima need two tables and
    pair-conditioned minima three; lazy builds omit the unused PF table for
    fresh minimum-only requests. Enumeration additionally reserves its frontier.
    Keep mask/output space and linear bookkeeping in every estimate. The
    three-table default matches the ABI-1 constructor's native cap exactly.
    """
    if type(n) is not int or not 0 < n <= 100000:
        raise ValueError("Expected 1–100000 nucleotides for native workspace sizing")
    if type(tables) is not int or not 1 <= tables <= 3:
        raise ValueError("Expected 1–3 connected-DP workspace tables")
    cells = ((n + 2) * (n + 1) // 2) * 21
    return cells * tables * 8 + n * n * 9 + (n + 1) * 256


@lru_cache(maxsize=1)
def _mach_api():
    # Public macOS host_statistics(HOST_VM_INFO=2), VM_STATISTICS_COUNT=15.
    # All counters here are natural_t (32-bit); only page counts are used.
    library = C.CDLL('/usr/lib/libSystem.B.dylib')
    library.mach_host_self.argtypes, library.mach_host_self.restype = [], C.c_uint32
    library.host_statistics.argtypes = [C.c_uint32, C.c_int, C.c_void_p, C.POINTER(C.c_uint32)]
    library.host_statistics.restype = C.c_int
    # Retain one host send right for this process, rather than obtaining a new
    # right at every complex evaluation.
    return library, library.mach_host_self()


def _mac_available() -> int:
    library, host = _mach_api()
    values, count = (C.c_uint32 * 15)(), C.c_uint32(15)
    # POSIX returns the page size as a Python integer and avoids a Mach
    # vm_size_t pointer (64-bit on both current macOS architectures).
    page_size = os.sysconf('SC_PAGESIZE')
    if (library.host_statistics(host, 2, C.byref(values), C.byref(count))
            or count.value < 3 or page_size <= 0):
        raise OSError("macOS memory statistics unavailable")
    # free_count already includes speculative pages. Purgeable pages overlap
    # existing categories; adding either count again would overstate headroom.
    return (values[0] + values[2]) * page_size


def _linux_available() -> int:
    data = {}
    for line in Path('/proc/meminfo').read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[0] in ('MemAvailable:', 'MemFree:'):
            data[fields[0][:-1]] = int(fields[1]) * 1024
    available = data.get('MemAvailable', data.get('MemFree'))
    if available is None:
        raise OSError("Linux memory statistics unavailable")
    # WSL reports the guest memory here. Containers may have a smaller cgroup
    # allowance; include both the mounted root and this process's membership.
    roots = [(Path('/sys/fs/cgroup'), 'memory.max', 'memory.current'),
             (Path('/sys/fs/cgroup/memory'), 'memory.limit_in_bytes', 'memory.usage_in_bytes')]
    try:
        for line in Path('/proc/self/cgroup').read_text().splitlines():
            _, controllers, relative = line.split(':', 2)
            path = PurePosixPath(relative)
            if '..' in path.parts:
                continue
            if not controllers:
                base, maximum, current = Path('/sys/fs/cgroup'), 'memory.max', 'memory.current'
            elif 'memory' in controllers.split(','):
                base, maximum, current = Path('/sys/fs/cgroup/memory'), 'memory.limit_in_bytes', 'memory.usage_in_bytes'
            else:
                continue
            member = base / str(path).lstrip('/')
            # A parent can impose the effective limit even when the leaf says
            # "max". Walk every containing cgroup up to the mounted root.
            while member != base:
                roots.append((member, maximum, current))
                member = member.parent
    except (OSError, ValueError):
        pass
    for root, maximum, current in dict.fromkeys(roots):
        try:
            limit = (root / maximum).read_text().strip()
            if limit != 'max':
                remaining = max(0, int(limit) - int((root / current).read_text().strip()))
                available = min(available, remaining)
        except (OSError, ValueError):
            continue
    return available


def available_memory_bytes() -> int | None:
    """Fresh OS observation without subprocesses or optional dependencies."""
    try:
        value = _mac_available() if sys.platform == 'darwin' else _linux_available() if sys.platform.startswith('linux') else None
        return value if type(value) is int and value >= 0 else None
    except (OSError, AttributeError, ValueError):
        return None


def safe_budget_bytes(available: int | None) -> int:
    if available is None:
        return _UNKNOWN_BUDGET
    return max(0, available - max(256 * _MIB, available // 4))


@contextmanager
def reserve_workspace(n: int, *, explicit_limit: int | None, output_bytes: int = 0,
                      python_bytes: int = 0, tables: int = 3,
                      enumerate_minima: bool = False, lazy_minimum: bool = False):
    """Admit one complete request before native/Python quadratic allocations.

    Explicit limits remain exact native caps. Automatic limits fit the required
    workspace and outputs rather than silently truncating the complex catalog.
    Host reservations count only the tables used by the observable, while the
    native cap still satisfies the constructor's three-table ABI-1 check.
    Python conversion allowance and concurrent reservations count against the
    host budget, in addition to the native workspace. No sequences or physical
    parameters change when admission fails.
    """
    global _RESERVED
    required = workspace_bytes(n)
    host_workspace = workspace_bytes(n, tables=tables)
    if type(enumerate_minima) is not bool:
        raise ValueError("enumerate_minima must be bool")
    if type(lazy_minimum) is not bool:
        raise ValueError("lazy_minimum must be bool")
    if explicit_limit is not None and max(required, output_bytes) > explicit_limit:
        raise MemoryError(f"Native workspace requires {required:,} bytes for {n} nucleotides; "
                          f"the explicit limit is {explicit_limit:,} bytes.")
    # The native estimate includes its largest fixed workspace. Add output
    # buffers and a small frontier allowance for bounded MFE enumeration.
    native_needed = required + output_bytes + 8 * _MIB
    native_limit = explicit_limit if explicit_limit is not None else native_needed
    reservation = host_workspace + output_bytes + 8 * _MIB + python_bytes
    if enumerate_minima:
        # connected_dp.cpp bounds its enumeration frontier with fixed_bytes =
        # memory_ - one DP table + bookkeeping. That bound charges an unused
        # n*n double output matrix: enumeration returns strings, not pairs.
        # Cover the entire permitted frontier (including a larger explicit cap),
        # subtract only that absent matrix, and add outputs plus a fresh margin
        # for parameter snapshots and temporary frontier growth. Two tables
        # alone would undercount the frontier, even for maximum_count=1.
        # New builds allocate the PF table only on its first use. The fresh
        # minimum-only handle therefore also omits that whole table, although
        # the native frontier guard still charges it in fixed_bytes.
        absent_inside = ((n + 2) * (n + 1) // 2) * 21 * 8 if lazy_minimum else 0
        reservation = max(reservation, native_limit - 8 * n * n - absent_inside
                          + output_bytes + 8 * _MIB + python_bytes)
    with _LOCK:
        available = available_memory_bytes()
        budget = safe_budget_bytes(available)
        if reservation > budget - _RESERVED:
            source = 'available system memory' if available is not None else 'unavailable system memory statistics (conservative budget)'
            raise MemoryError(
                f"This complex contains {n} nucleotides and requires approximately "
                f"{reservation / 1024**3:.2f} GiB of workspace/output memory; "
                f"{max(0, budget - _RESERVED) / 1024**3:.2f} GiB is currently available "
                f"after reserving host headroom ({source}). Close other calculations "
                "or use a machine with more available memory; the requested complex catalog is unchanged.")
        _RESERVED += reservation
    try:
        yield native_limit
    finally:
        with _LOCK:
            _RESERVED -= reservation
