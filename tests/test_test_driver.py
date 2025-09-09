#!/usr/bin/python

import os
import subprocess

import kim_edn
import numpy as np
import numpy.typing as npt
from ase.atoms import Atoms
from ase.build import bulk
from ase.calculators.lj import LennardJones

from kim_tools import (
    KIMTestDriver,
    SingleCrystalTestDriver,
    cartesian_rotation_is_in_point_group,
    detect_unique_crystal_structures,
    get_deduplicated_property_instances,
)
from kim_tools.test_driver.core import _get_optional_source_value


class TestInitSingleCrystalTestDriver(SingleCrystalTestDriver):
    def _calculate(self, **kwargs) -> None:
        pass


class TestInitKIMTestDriver(KIMTestDriver):
    def _calculate(self, **kwargs) -> None:
        pass


class TestIsolatedEnergyDriver(KIMTestDriver):
    def _calculate(self, species):
        """
        Example calculate method for testing isolated energy getter
        """
        assert np.isclose(self.get_isolated_energy_per_atom(species), 0.0)


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


class TestStructureDetectionTestDriver(SingleCrystalTestDriver):
    def _calculate(self, deform_matrix: npt.NDArray = np.eye(3), **kwargs):
        """
        strain the crystal and write a crystal-structure-npt.
        For testing various crystal detection things
        """
        atoms = self._get_atoms()
        original_cell = atoms.cell
        new_cell = original_cell @ deform_matrix
        atoms.set_cell(new_cell, scale_atoms=True)
        self._update_nominal_parameter_values(atoms)
        self._add_property_instance_and_common_crystal_genome_keys(
            "crystal-structure-npt"
        )


class FileWritingTestDriver(KIMTestDriver):
    def _calculate(self):
        """
        Mock calculate for file writing testing
        """
        self._add_property_instance("file-prop")
        with open("foo.txt", "w") as f:
            f.write("foo")
        self._add_file_to_current_property_instance("textfile", "foo.txt")


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


def test_structure_detection():
    test = TestStructureDetectionTestDriver(LennardJones())
    atoms = bulk("Mg")
    hcp_prototype = "A_hP2_194_c"
    hcp_library_prototype = "A_hP2_194_c-001"
    hcp_shortname = ["Hexagonal Close Packed (Mg, $A3$, hcp) Structure"]
    stretch = np.diag([1, 1, 10])
    for deform, prototype, library_prototype, shortname in zip(
        [np.eye(3), stretch],
        [hcp_prototype, hcp_prototype],
        [hcp_library_prototype, None],
        [hcp_shortname, None],
    ):
        property_instance = test(atoms, deform_matrix=deform)[0]
        assert property_instance["prototype-label"]["source-value"] == prototype
        assert (
            _get_optional_source_value(property_instance, "library-prototype-label")
            == library_prototype
        )
        assert _get_optional_source_value(property_instance, "short-name") == shortname


def test_get_isolated_energy_per_atom():
    for model in [
        LennardJones(),
        "LennardJones612_UniversalShifted__MO_959249795837_003",
        "Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu",
    ]:
        td = TestIsolatedEnergyDriver(model)
        for species in ["Pb", "Au"]:
            td(species=species)


def test_init_rotation():
    td = TestInitSingleCrystalTestDriver(LennardJones())
    atoms = bulk("Mg", crystalstructure="bct", a=1, c=3.14)
    atoms.rotate(60, "z", rotate_cell=True)
    # We rotated from standard orientation by 60 deg
    # ccw. So when the atoms get rebuilt, the rotation
    # we will do to standardize is 60 deg cw.
    ref_rotation = np.array(
        [  # 60 deg clockwise
            [0.5, 0.866025, 0.0],
            [-0.866025, 0.5, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    td(atoms)
    input_rotation = td.get_input_rotation()
    assert cartesian_rotation_is_in_point_group(
        ref_rotation.T @ input_rotation,
        139,
        atoms.cell,
    )


def test_file_writing():
    td = FileWritingTestDriver(LennardJones())
    td()
    td.write_property_instances_to_file()
    assert os.path.isfile("output/foo-1.txt")
    # add a file into output to make sure it's not being moved
    with open("output/bar", "w") as f:
        f.write("bar")
    td.write_property_instances_to_file("foo.edn")
    # will raise a filenotfound error if missing
    os.remove("foo.edn")
    os.remove("foo-1.txt")
    os.remove("output/bar")


def test_atom_style():
    td = TestInitKIMTestDriver("LennardJones612_UniversalShifted__MO_959249795837_003")
    assert td._get_supported_lammps_atom_style() == "atomic"
    charge_sm = "Sim_LAMMPS_ReaxFF_AnGoddard_2015_BC__SM_389039364091_000"
    subprocess.run(
        f"kim-api-collections-management install --force CWD {charge_sm}",
        shell=True,
        check=True,
    )
    td = TestInitKIMTestDriver(charge_sm)
    assert td._get_supported_lammps_atom_style() == "charge"
    subprocess.run(
        f"kim-api-collections-management remove --force {charge_sm}",
        shell=True,
        check=True,
    )


if __name__ == "__main__":
    test_atom_style()
