import subprocess
from warnings import warn

import setuptools

try:
    subprocess.check_output(["aflow", "--proto=A_cF4_225_a"])
except Exception:
    message = (
        "aflow executable not found in PATH. "
        "You will not be able to run any Crystal Genome tests."
    )
    lines = "=" * 89
    warn(message)
    print()
    print(lines)
    print(message)
    print(lines)
    print()

try:
    subprocess.check_output(["units", "--help"])
except Exception:
    message = (
        "GNU `units` executable not found in PATH. "
        "Unit conversions will not be available."
    )
    lines = "=" * 89
    warn(message)
    print()
    print(lines)
    print(message)
    print(lines)
    print()


setuptools.setup()
