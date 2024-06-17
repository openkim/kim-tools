import setuptools
import subprocess
from warnings import warn

with open("README.md", "r") as fh:
    long_description = fh.read()

try:
    subprocess.check_output('aflow --proto=A_cF4_225_a')
except Exception:
    warn("aflow executable not found in PATH. You will not be able to run any Crystal Genome tests.")

setuptools.setup(
    name="kim-test-utils",
    version="0.1.0",
    description=(
        "Helper routines for writing KIM Tests"
    ),
    author=["ilia Nikiforov, Eric Fuemmler, Ellad Tadmor"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    license="CDDL",
    install_requires=["ase >= 3.19.0b1", "kim-property >= 2.6.2", "kim-edn >= 1.4.1", "kimpy >= 2.1.1", "spglib >= 2.1.0"],
    classifiers=[
        "Development Status :: 4 - Beta"
        "License :: OSI Approved :: Common Development and Distribution License 1.0 (CDDL-1.0)",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    include_package_data=True
)
