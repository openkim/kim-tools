#!/usr/bin/python

import os

import kim_edn
from ase.atoms import Atoms
from ase.calculators.lj import LennardJones

from kim_tools import (
    KIMTestDriver,
    detect_unique_crystal_structures,
    get_deduplicated_property_instances,
)


class TestTestDriver(KIMTestDriver):
    def _calculate(self, property_name, species):
        """
        example calculate method

        Args:
            property_name: for testing ability to find properties at different paths.
            !!! AN ACTUAL TEST DRIVER SHOULD NOT HAVE AN ARGUMENT SUCH AS THIS !!!
        """
        atoms = Atoms([species], [[0, 0, 0]])
        self._add_property_instance(property_name, "This is an example disclaimer.")
        self._add_key_to_current_property_instance(
            "species", atoms.get_chemical_symbols()[0]
        )
        self._add_key_to_current_property_instance(
            "mass", atoms.get_masses()[0], "amu", {"source-std-uncert-value": 1}
        )


def test_kimtest(monkeypatch):
    test = TestTestDriver(LennardJones())
    testing_property_names = [
        "atomic-mass",  # already in kim-properties
        "atomic-mass0",  # found in $PWD/local-props
        "atomic-mass0",  # check that repeat works fine
        # check that full id works as well, found in $PWD/local-props
        "tag:brunnels@noreply.openkim.org,2016-05-11:property/atomic-mass1",
        "atomic-mass2",  # found in $PWD/local-props/atomic-mass2
        "atomic-mass3",  # found in $PWD/mock-test-drivers-dir/mock-td/local_props,
        # tested using the monkeypatch below
    ]

    monkeypatch.setenv(
        "KIM_PROPERTY_PATH",
        os.path.join(os.getcwd(), "mock-test-drivers-dir/*/local-props")
        + ":"
        + os.path.join(os.getcwd(), "mock-test-drivers-dir/*/local_props"),
    )

    for prop_name in testing_property_names:
        test(property_name=prop_name, species="Ar")

    assert len(test.property_instances) == 6
    test.write_property_instances_to_file()


def test_detect_unique_crystal_structures():
    reference_structure = kim_edn.load("structures/OSi.edn")
    test_structure = kim_edn.load("structures/OSi_twin.edn")
    assert (
        len(
            detect_unique_crystal_structures(
                [
                    reference_structure,
                    reference_structure,
                    test_structure,
                    test_structure,
                    test_structure,
                    test_structure,
                ],
                allow_rotation=True,
            )
        )
        == 1
    )
    assert (
        len(
            detect_unique_crystal_structures(
                [
                    reference_structure,
                    reference_structure,
                    test_structure,
                    test_structure,
                    test_structure,
                    test_structure,
                ],
                allow_rotation=False,
            )
        )
        == 2
    )


def test_get_deduplicated_property_instances():
    property_instances = kim_edn.load("structures/results.edn")
    fully_deduplicated = get_deduplicated_property_instances(property_instances)
    assert len(fully_deduplicated) == 6
    inst_with_1_source = 0
    inst_with_2_source = 0
    for property_instance in fully_deduplicated:
        n_inst = len(
            property_instance["crystal-genome-source-structure-id"]["source-value"][0]
        )
        if n_inst == 1:
            inst_with_1_source += 1
        elif n_inst == 2:
            inst_with_2_source += 1
        else:
            assert False
    assert inst_with_1_source == 3
    assert inst_with_2_source == 3
    partially_deduplicated = get_deduplicated_property_instances(
        property_instances, ["mass-density-crystal-npt"]
    )
    assert len(partially_deduplicated) == 8
    inst_with_1_source = 0
    inst_with_2_source = 0
    for property_instance in partially_deduplicated:
        n_inst = len(
            property_instance["crystal-genome-source-structure-id"]["source-value"][0]
        )
        if n_inst == 1:
            inst_with_1_source += 1
        elif n_inst == 2:
            inst_with_2_source += 1
        else:
            assert False
    assert inst_with_1_source == 7
    assert inst_with_2_source == 1


if __name__ == "__main__":
    test_get_deduplicated_property_instances()
