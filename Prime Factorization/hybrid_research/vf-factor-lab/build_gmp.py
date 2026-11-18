"""Optional Python convenience wrapper; Makefile.gmp builds without Python."""
from pathlib import Path
import os
import shlex
import subprocess

ROOT = Path(__file__).resolve().parent
out = ROOT / "build" / "vf_factor_gmp"
out.parent.mkdir(exist_ok=True)
command = (shlex.split(os.environ.get("CXX", "g++"))
           + shlex.split(os.environ.get("CPPFLAGS", ""))
           + shlex.split(os.environ.get("CXXFLAGS", "-O3 -std=c++17 -Wall -Wextra"))
           + [str(ROOT/"cpp/gmp_factor.cpp"), str(ROOT/"vf_factor/fusion_kernel.cpp")]
           + shlex.split(os.environ.get("LDFLAGS", ""))
           + ["-lgmpxx", "-lgmp", "-o", str(out)])
subprocess.run(command, check=True, cwd=ROOT)
subprocess.run([str(out), "--metadata"], check=True, cwd=ROOT)
print(out)
