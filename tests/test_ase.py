import numpy as np
from ase.calculators.lj import LennardJones

from kim_tools import get_isolated_energy_per_atom

# from lj_fail_no_neighbors import LennardJonesFailNoNeighbors


def test_get_isolated_energy_per_atom():
    assert np.isclose(
        get_isolated_energy_per_atom(
            "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "H"
        ),
        0,
    )
    assert np.isclose(
        get_isolated_energy_per_atom(LennardJones(), "H"),
        0,
    )
    # assert np.isclose(
    #     get_isolated_energy_per_atom(LennardJonesFailNoNeighbors(), "H"),
    #    0,
    # )
