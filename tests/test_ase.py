import numpy as np
from ase.calculators.kim import KIM
from ase.calculators.lj import LennardJones
from lj_fail_no_neighbors import LennardJonesFailNoNeighbors

from kim_tools import get_isolated_energy_per_atom

MO_NAME = "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
SM_NAME = "Sim_LAMMPS_ADP_StarikovGordeevLysogorskiy_2020_SiAuAl__SM_113843830602_000"


def test_get_isolated_energy_per_atom():
    for model in [
        LennardJones(),
        MO_NAME,
        SM_NAME,
        KIM(SM_NAME),  # This creates a LAMMPSLib object
        LennardJonesFailNoNeighbors(),  # This intentionally crashes for isolated atoms
    ]:
        for species in ["Au", "Al"]:
            assert np.isclose(
                get_isolated_energy_per_atom(model=model, symbol=species),
                0,
            )
