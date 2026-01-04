# kim-tools

[![Testing](https://github.com/openkim/kim-tools/actions/workflows/testing.yml/badge.svg)](https://github.com/openkim/kim-tools/actions/workflows/testing.yml)
[![docs](https://app.readthedocs.org/projects/kim-tools/badge/?version=latest)](https://kim-tools.readthedocs.io/en/latest/)
[![PyPI](https://img.shields.io/pypi/v/kim-tools.svg)](https://pypi.org/project/kim-tools/)
[![codecov](https://codecov.io/gh/openkim/kim-tools/graph/badge.svg?token=G57VDZYY0F)](https://codecov.io/gh/openkim/kim-tools)

KIMTestDriver and SingleCrystalTestDriver classes for creating OpenKIM Test Drivers, and helper routines for writing
KIM Tests and Verification Checks. Documentation at https://kim-tools.readthedocs.io.

## Contributing Guide (Under Construction)

All contributed functions, classes and methods should be documented with Google style docstrings (https://google.github.io/styleguide/pyguide.html#383-functions-and-methods) and should have type hints for all arguments and return values (https://docs.python.org/3/library/typing.html). The docstrings are automatically rendered into the API documentation. To check that this is working correctly, install the Python packages in `docs/requirements-docs.txt` and run `make html` in `docs/`. There are some warnings, but there should be no errors, and the rendered docstring in the appropriate module html (e.g., `docs/build/html/kim_tools.ase.html`) should have a properly rendered section for the functions, classes and methods you wrote. See `kim_tools/test_driver/core.py` for a variety of examples of docstrings.

The code has a simple test suite using [pytest](https://docs.pytest.org/en/stable/). See the various files named `test_*` in `tests/` for examples and add tests for any code you write.

To confirm that your code will pass the linter and style checks, install pre-commit as in this guide: https://pre-commit.com/#quick-start. The checks will run every time you make a commit.

`kim_tools/ase` contains utility functions that can be used for interatomic calculations. It is inherited from merging https://github.com/openkim/kim-python-utils into this repo. Previously, it worked only with KIM Models. All new functions added to this module should support being passed both a KIM Model as a string, and an ASE `Calculator` object. See `kim-tools/ase/core.py::get_isolated_energy_per_atom` for an example. If you are working on or with existing functions in this module, try to upgrade it to this functionality.
