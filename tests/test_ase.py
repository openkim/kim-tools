import numpy as np
from ase.calculators.lj import LennardJones

from kim_tools import get_isolated_energy_per_atom

# from lj_fail_no_neighbors import LennardJonesFailNoNeighbors


def test_get_isolated_energy_per_atom():
    for model in [
        LennardJones(),
        "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003",
        "Sim_LAMMPS_Table_GrogerVitekDlouhy_2020_CoCrFeMnNi__SM_786004631953_001",
    ]:
        for species in ["Fe", "Cr"]:
            assert np.isclose(
                get_isolated_energy_per_atom(model=model, symbol=species),
                0,
            )
    """
    assert np.isclose(
         get_isolated_energy_per_atom(LennardJonesFailNoNeighbors(), "H"),
        0,
     )
     """
