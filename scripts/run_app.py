from __future__ import annotations

import sys


if __name__ == "__main__":
    from aptaswitch_core.nupack_setup import handle_nupack_probe

    if handle_nupack_probe():
        raise SystemExit(0)
    if len(sys.argv) == 3 and sys.argv[1] == "--smoke-test":
        from aptaswitch_studio.smoke_check import run_smoke_check

        raise SystemExit(run_smoke_check(sys.argv[2]))
    from aptaswitch_studio.__main__ import main

    raise SystemExit(main())
