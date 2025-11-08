import json
import subprocess
from os import PathLike
from typing import Dict, List, Optional

import pytest

from kim_tools import detect_unique_crystal_structures, query_crystal_structures

# TEST_CASES = [572, 1656]  # Test just on triclinic and cubic, extremes of symmetry
TEST_CASES = [572, 365, 1729, 1194, 1473, 166, 1205, 1357, 915, 212, 641, 22]
MATERIALS_FILE = "test_structures.json"
QUERY_DUMP = "output/query_result.json"


def get_test_crystal_structures(
    materials_file: PathLike = MATERIALS_FILE,
    test_cases: Optional[List[int]] = TEST_CASES,
    deduplicate: bool = False,
) -> List[Dict]:
    """
    Query OpenKIM reference data for materials (prototype label + species) from
    `materials_file`. Deduplicate the structures if asked, otherwise just return the
    first index of each (you don't want to be testing 25 identical copies of FCC
    Aluminum)

    Returns a list of `crystal-structure-npt` property instances
    """
    with open(materials_file) as f:
        query_inputs = json.load(f)

    if test_cases is not None:
        query_inputs = [query_inputs[i] for i in test_cases]

    test_crystal_structures = []
    for query_input in query_inputs:
        query_result = query_crystal_structures(**query_input)
        if deduplicate:
            indices_to_test = detect_unique_crystal_structures(query_result)
        else:
            indices_to_test = [0]
        test_crystal_structures += [query_result[i] for i in indices_to_test]

    with open(QUERY_DUMP, "w") as f:
        json.dump(test_crystal_structures, f)

    return test_crystal_structures


@pytest.fixture(scope="session")
def input_crystal_structures():
    return get_test_crystal_structures()


@pytest.fixture(scope="session", autouse=True)
def install_models():
    try:
        subprocess.run(
            [
                "kim-api-collections-management",
                "install",
                "CWD",
                "Sim_LAMMPS_ReaxFF_AnGoddard_2015_BC__SM_389039364091_000",
            ],
            check=True,
            capture_output=True,
            encoding="utf-8",
        )
    except subprocess.CalledProcessError as e:
        if ("already installed" in e.output) or (
            "does not appear to contain CMakeLists.txt." in e.output
        ):
            pass
        else:
            raise e
    try:
        subprocess.run(
            [
                "kim-api-collections-management",
                "install",
                "CWD",
                "LennardJones612_UniversalShifted__MO_959249795837_003",
            ],
            check=True,
            capture_output=True,
            encoding="utf-8",
        )
    except subprocess.CalledProcessError as e:
        if ("already installed" in e.output) or (
            "does not appear to contain CMakeLists.txt." in e.output
        ):
            pass
        else:
            raise e
