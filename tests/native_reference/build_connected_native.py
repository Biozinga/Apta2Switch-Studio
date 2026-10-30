"""Build the optional synthetic C++ grammar library without installing anything."""
import argparse
from pathlib import Path
import platform
import subprocess
from tempfile import TemporaryDirectory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compiler", default="clang++")
    parser.add_argument("--sanitize", action="store_true", help="Enable address/undefined behavior sanitizers")
    parser.add_argument("--check-factory", action="store_true", help="Compile and run standalone C++ provider checks")
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[2] / "src/aptaswitch_core/thermodynamics/native/connected_dp.cpp"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    command = [args.compiler, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror", "-fPIC",
               "-dynamiclib" if platform.system() == "Darwin" else "-shared"]
    if args.sanitize:
        command += ["-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
    command += [str(source), "-o", str(args.output)]
    subprocess.run(command, check=True)
    if args.check_factory:
        with TemporaryDirectory(prefix="aptaswitch-native-factory-") as temporary:
            executable = Path(temporary) / "factory-check"
            check = [args.compiler, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror"]
            if args.sanitize:
                check += ["-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
            check += ["-I", str(source.parent), str(source),
                      str(Path(__file__).with_name("test_connected_native_factory.cpp")), "-o", str(executable)]
            subprocess.run(check, check=True)
            subprocess.run([str(executable)], check=True, timeout=30)
    print(args.output.resolve())


if __name__ == "__main__":
    main()
