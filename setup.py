import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

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
    install_requires=["ase >= 3.19.0b1", "kim-property >= 2.5.8"],
    classifiers=[
        "Development Status :: 4 - Beta"
        "License :: OSI Approved :: Common Development and Distribution License 1.0 (CDDL-1.0)",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    include_package_data=True
)
