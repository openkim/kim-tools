import numpy as np
from ase.calculators.kim import KIM
from ase.calculators.lj import LennardJones
from lj_fail_no_neighbors import LennardJonesFailNoNeighbors

from kim_tools import find_equilibrium_config_FCC, get_isolated_energy_per_atom

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


# end-to-end test
def test_find_equilibrium_configuration_FCC():
    model = "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
    species = ["Au", "Al"]
    result = find_equilibrium_config_FCC(model, species)

    assert np.isclose(result["mono_species_equilibrium_alats"]["Al"], 3.3281, atol=1e-4)
    assert np.isclose(result["mono_species_equilibrium_alats"]["Au"], 3.7407, atol=1e-4)
    assert np.isclose(result["equilibrium_alat"], 3.5630, atol=1e-4)
    assert np.isclose(result["approx_mixed_equilibrium_alat"], 3.5344, atol=1e-4)
    assert result["ncells_per_side"] == 1


# end-to-end test (will generate an FCC config larger than default)
def test_find_equilibrium_configuration_FCC_n_cells_per_side():
    model = "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
    species = ["H", "C", "N", "O", "S"]
    result = find_equilibrium_config_FCC(model, species)

    assert np.isclose(result["mono_species_equilibrium_alats"]["H"], -1, atol=1e-4)
    assert np.isclose(result["mono_species_equilibrium_alats"]["C"], 2.0904, atol=1e-4)
    assert np.isclose(result["mono_species_equilibrium_alats"]["N"], 1.9528, atol=1e-4)
    assert np.isclose(result["mono_species_equilibrium_alats"]["O"], 1.8153, atol=1e-4)
    assert np.isclose(result["mono_species_equilibrium_alats"]["S"], 2.8880, atol=1e-4)
    assert np.isclose(result["equilibrium_alat"], 2.4201, atol=1e-4)
    assert np.isclose(result["approx_mixed_equilibrium_alat"], 2.1866, atol=1e-4)
    assert result["ncells_per_side"] == 2


# expected to fail due to coarse-scanning timeout
def test_find_equilibrium_configuration_FCC_coarseTimeout():
    model = "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
    species = ["Au", "Al"]
    result = find_equilibrium_config_FCC(model, species, coarse_timeout=1e-16)

    assert result["mono_species_equilibrium_alats"]["Al"] == -1.0
    assert result["mono_species_equilibrium_alats"]["Au"] == -1.0
    assert result["equilibrium_alat"] == -1.0
    assert np.isclose(result["approx_mixed_equilibrium_alat"], -1.0, atol=1e-4)
    assert result["mono_species_results"][0]["reference_config"]["ncells_per_side"] == 1
    assert result["mono_species_results"][1]["reference_config"]["ncells_per_side"] == 1
    assert result["mixed_species_result"]["reference_config"]["ncells_per_side"] == 1


# expected to fail due to nelder-mead timeout
def test_find_equilibrium_configuration_FCC_nelderMeadTimeout():
    model = "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
    species = ["Au", "Al"]
    result = find_equilibrium_config_FCC(model, species, nelder_mead_timeout=1e-16)

    assert result["mono_species_equilibrium_alats"]["Al"] == -1.0
    assert result["mono_species_equilibrium_alats"]["Au"] == -1.0
    assert result["equilibrium_alat"] == -1.0
    assert np.isclose(result["approx_mixed_equilibrium_alat"], -1.0, atol=1e-4)
    assert result["mono_species_results"][0]["reference_config"]["ncells_per_side"] == 1
    assert result["mono_species_results"][1]["reference_config"]["ncells_per_side"] == 1
    assert result["mixed_species_result"]["reference_config"]["ncells_per_side"] == 1
