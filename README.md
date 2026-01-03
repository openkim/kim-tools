# kim-tools

[![Testing](https://github.com/openkim/kim-tools/actions/workflows/testing.yml/badge.svg)](https://github.com/openkim/kim-tools/actions/workflows/testing.yml)
[![docs](https://app.readthedocs.org/projects/kim-tools/badge/?version=latest)](https://kim-tools.readthedocs.io/en/latest/)
[![PyPI](https://img.shields.io/pypi/v/kim-tools.svg)](https://pypi.org/project/kim-tools/)
[![codecov](https://codecov.io/gh/openkim/kim-tools/graph/badge.svg?token=G57VDZYY0F)](https://codecov.io/gh/openkim/kim-tools)

KIMTestDriver and SingleCrystalTestDriver classes for creating OpenKIM Test Drivers, and helper routines for writing
KIM Tests and Verification Checks. Documentation at https://kim-tools.readthedocs.io.

## Contributing Guide (Under Construction)

All contributed functions and methods should be documented with Google style docstrings (https://google.github.io/styleguide/pyguide.html#383-functions-and-methods) and should have type hints for all arguments and return values (https://docs.python.org/3/library/typing.html).

The code has a simple test suite using [pytest](https://docs.pytest.org/en/stable/). See the various files named `test_*` in `tests/` for examples and add tests for any code you write.

When contributing to `kim-tools/ase`, all new functions should support being passed both a KIM Model as a string, and an ASE `Calculator` object. See `kim-tools/ase/core.py::get_isolated_energy_per_atom` for an example. If you are working on an existing function in this module, try to upgrade it to this functionality.
