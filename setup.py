import setuptools
import subprocess
from warnings import warn

try:
    subprocess.check_output(['aflow','--proto=A_cF4_225_a'])
except Exception:
    message = "aflow executable not found in PATH. You will not be able to run any Crystal Genome tests."
    lines =   "========================================================================================="
    warn(message)
    print()
    print(lines)
    print(message)
    print(lines)
    print()

try:
    subprocess.check_output(['units','--help'])
except Exception:
    raise RuntimeError("units executable not found. It is required.")

setuptools.setup()
