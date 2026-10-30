"""Bundled isolated Python/oxpy process; speaks the existing worker protocol."""
from pathlib import Path
import runpy
import sys


if __name__ == "__main__":
    if sys.argv[1:2] == ["-c"]:
        exec(compile(sys.argv[2], "<oxpy-runtime-probe>", "exec"), {"__name__": "__main__"})
    else:
        script = Path(sys.argv[1]).resolve()
        if script.name not in {"oxdna_ffs_worker.py", "oxdna_kinetics_worker.py", "oxdna_native_ffs.py", "oxdna_guided_worker.py", "oxdna_free_md_worker.py"}:
            raise SystemExit("Unsupported bundled oxDNA worker")
        sys.argv = sys.argv[1:]
        sys.path.insert(0, str(script.parent))
        runpy.run_path(str(script), run_name="__main__")
