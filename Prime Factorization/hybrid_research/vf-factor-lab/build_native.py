#!/usr/bin/env python3
"""Build optional native64 fairness baselines with a GCC/Clang-style compiler.

Usage: python3 build_native.py [--compiler g++]
Requires a C++17 compiler supporting unsigned __int128 (GCC or Clang).
No external C++ libraries are required. Native binaries are platform-specific.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import platform
import shlex
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default=os.environ.get("CXX", "g++"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    source = root / "vf_factor" / "native_kernel.cpp"
    system = platform.system()
    suffix = ".dylib" if system == "Darwin" else ".dll" if system == "Windows" else ".so"
    output = source.with_suffix(suffix)
    linker = "-dynamiclib" if system == "Darwin" else "-shared"
    command = shlex.split(args.compiler) + [
        "-O3", "-std=c++17", "-Wall", "-Wextra", "-fPIC", linker,
        str(source), "-o", str(output),
    ]
    subprocess.run(command, check=True)
    print(output)
    from vf_factor.native import backend_info
    print(backend_info())


if __name__ == "__main__":
    main()
