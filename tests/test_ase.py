import numpy as np
from ase.calculators.kim import KIM
from ase.calculators.lj import LennardJones
from lj_fail_no_neighbors import LennardJonesFailNoNeighbors

from kim_tools import (
    fcc_atoms_in_supercell,
    find_equilibrium_config_FCC,
    generate_fcc_compute_energy,
    get_isolated_energy_per_atom,
)

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


def test_fcc_atoms_in_supercell():
    assert fcc_atoms_in_supercell(1) == 4
    assert fcc_atoms_in_supercell(2) == 32
    assert fcc_atoms_in_supercell(3) == 108


def test_generate_fcc_compute_energy():

    # test with model name
    model = MO_NAME
    species = ["Au", "Al"]
    alat = 3.5
    seed = 13
    energy, ncells = generate_fcc_compute_energy(model, species, alat, seed)
    assert np.isclose(energy, -614.5593, atol=1e-4)
    assert np.isclose(ncells, 2)

    # test with KIM calculator object
    # model = KIM(SM_NAME)
    model = SM_NAME
    species = ["Au", "Al"]
    alat = 3.5
    seed = 13
    energy, ncells = generate_fcc_compute_energy(model, species, alat, seed)
    assert np.isclose(energy, 14.3366, atol=1e-4)
    assert np.isclose(ncells, 2)


# end-to-end test
def test_find_working_configuration_FCC():
    model = "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
    species = ["Au", "Al"]
    result = find_equilibrium_config_FCC(model, species)

    assert np.isclose(result["mono_species_equilibrium_alats"]["Al"], 3.3281, atol=1e-4)
    assert np.isclose(result["mono_species_equilibrium_alats"]["Au"], 3.7407, atol=1e-4)
    assert np.isclose(result["equilibrium_alat"], 3.5768, atol=1e-4)
    assert result["ncells_per_side"] == 2
