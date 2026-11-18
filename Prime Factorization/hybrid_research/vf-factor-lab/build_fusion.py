#!/usr/bin/env python3
"""Build the resumable native wheel cursor used by the unified factorizer."""
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
    source = Path(__file__).resolve().parent / "vf_factor" / "fusion_kernel.cpp"
    system = platform.system()
    suffix = ".dylib" if system == "Darwin" else ".dll" if system == "Windows" else ".so"
    output = source.with_suffix(suffix)
    linker = "-dynamiclib" if system == "Darwin" else "-shared"
    subprocess.run(shlex.split(args.compiler) + [
        "-O3", "-std=c++17", "-Wall", "-Wextra", "-fPIC", linker,
        str(source), "-o", str(output),
    ], check=True)
    print(output)
    from vf_factor.fusion_native import backend_info
    print(backend_info())


if __name__ == "__main__":
    main()
