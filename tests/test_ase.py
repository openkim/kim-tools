import numpy as np
from ase.calculators.kim import KIM
from ase.calculators.lj import LennardJones
from lj_fail_no_neighbors import LennardJonesFailNoNeighbors

from kim_tools import (
    fcc_atoms_in_supercell,
    filter_good_alat,
    find_working_configuration_FCC,
    generate_fcc_compute_energy,
    get_isolated_energy_per_atom,
    local_edge_detection,
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
    assert np.isclose(energy, -288.4856, atol=1e-4)
    assert np.isclose(ncells, 2)

    # test with KIM calculator object
    model = KIM(SM_NAME)
    species = ["Au", "Al"]
    alat = 3.5
    seed = 13
    energy, ncells = generate_fcc_compute_energy(model, species, alat, seed)
    assert np.isclose(energy, -70.2935, atol=1e-4)
    assert np.isclose(ncells, 2)


def test_local_edge_detection():
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y = [1, 2, 1, 2, 1, 50, 1, 2, 1, 2]
    leds = local_edge_detection(x, y)
    expected_leds = [10.666, -42.666, 82.666, -82.666, 42.666]
    assert np.allclose(leds, expected_leds, atol=1e-3)

    y = np.sin(np.array(x))
    leds = local_edge_detection(x, y)
    expected_leds = [-0.1265, -0.0284, 0.0957, 0.1319, 0.0468]
    assert np.allclose(leds, expected_leds, atol=1e-3)


def test_filter_good_alat():
    alats = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    energies = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
    ncells = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    leds = [0.1, 0.5, 1.5, 0.8, 0.3, 2.0, 0.4, 0.2, 1.2, 0.6]
    min_cutoff = 0.5
    energy_bound = [5e-2, 5e2]
    led_tol = 1.0
    result = filter_good_alat(
        alats, energies, ncells, leds, min_cutoff, energy_bound, led_tol
    )

    # check that the return type is dict and contains expected keys
    assert isinstance(result, dict)
    assert len(result) == 3
    assert "good_alat" in result
    assert "min_led" in result
    assert "good_ncells" in result


# end-to-end test
def test_find_working_configuration_FCC():
    model = "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
    species = ["Au", "Al"]
    energy_bound = [5e-2, 5e2]
    led_tol = 1.0
    seed = 21
    result = find_working_configuration_FCC(model, species, energy_bound, led_tol, seed)

    assert np.isclose(result["good_alat"], 5.2071, atol=1e-4)
    assert result["good_ncells"] == 2
    assert np.isclose(result["min_led"], 3.6274390898446044e-11, atol=1e-15)
