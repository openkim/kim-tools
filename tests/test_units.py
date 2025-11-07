import numpy as np

from kim_tools.kimunits import convert_list, convert_units


def test_units() -> None:
    one_kilogram = convert_units(1000, "g")
    assert np.isclose(one_kilogram[0], 1)
    assert one_kilogram[1] == "kg"
    assert np.isclose(convert_units(0, "tempC", "tempF", suppress_unit=True), 32)
    assert np.allclose(convert_list([1000, 1000], "g")[0], [1, 1])
