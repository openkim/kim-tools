#!/usr/bin/python

import glob
import os
import shutil
import tarfile
from tempfile import TemporaryDirectory

import kim_edn
import numpy as np
import numpy.typing as npt
from ase.atoms import Atoms
from ase.build import bulk
from ase.calculators.lj import LennardJones

from kim_tools import (
    KIMTestDriver,
    KIMTestDriverError,
    SingleCrystalTestDriver,
    cartesian_rotation_is_in_point_group,
    detect_unique_crystal_structures,
    get_deduplicated_property_instances,
)
from kim_tools.test_driver.core import TOKENPATH, _get_optional_source_value


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
    def _calculate(self, write_aux_files: bool = True):
        """
        Mock calculate for file writing testing. Will result in the following
        files being created under the output directory and listed in the property
        instance (n is the number of times _calculate has been previously called):

        output/foo_prop-{3*n+1}.txt
        output/bar/bar_prop-{3*n+2}.txt
        output/baz_prop-{3*n+3}

        And an archive named aux_files.{n}.txz containing the following directory
        tree:
        aux_files.{n}/foo
        aux_files.{n}/bar/bar
        aux_files.{n}/baz/baz

        Args:
            write_aux_files:
                Whether to write the above-mentioned aux files
        """
        self._add_property_instance("file-prop")
        with open("output/foo_prop.txt", "w") as f:
            f.write("foo")
        # Because it's under output, should leave it untouched except rename
        # to foo_prop-{n+1}.txt (where n is the number of times _calculate has been
        # previously called)
        self._add_file_to_current_property_instance("textfile", "output/foo_prop.txt")

        self._add_property_instance("file-prop")
        os.makedirs("bar", exist_ok=True)
        with open("bar/bar_prop.txt", "w") as f:
            f.write("bar")
        # Not under output, but under CWD, so should preserve relative path and
        # write to "output/bar/bar_prop-{n+1}.txt (where n is the number of times
        # _calculate has been previously called)
        self._add_file_to_current_property_instance("textfile", "bar/bar_prop.txt")
        os.rmdir("bar")

        self._add_property_instance("file-prop")
        with TemporaryDirectory() as d:
            tempfilepath = os.path.join(d, "baz_prop")
            with open(tempfilepath, "w") as f:
                f.write("baz")
            # Random path in system not under CWD, so should just go to
            # "output/baz_prop-{n+1}.txt (where n is the number of times
            # _calculate has been previously called)
            self._add_file_to_current_property_instance("textfile", tempfilepath)

        if write_aux_files:
            # Now, write some auxiliary files
            with open("output/foo", "w") as f:
                f.write("foo")
            with open("output/bar/bar", "w") as f:
                f.write("bar")
            os.makedirs("output/baz", exist_ok=True)
            with open("output/baz/baz", "w") as f:
                f.write("baz")


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
        "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003",
        "Sim_LAMMPS_ADP_StarikovGordeevLysogorskiy_2020_SiAuAl__SM_113843830602_000",
    ]:
        td = TestIsolatedEnergyDriver(model)
        for species in ["Al", "Au"]:
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
    """
    The n-th call (zero-based) to FileWritingTestDriver should
    result in the following directory structure:

    output/foo_prop-{n+1}.txt
    output/bar/bar_prop-{n+2}.txt
    output/baz_prop-{n+3}

    And an archive named output/aux_files.{n}.txz containing the following directory
    tree:
    aux_files.{n}/foo
    aux_files.{n}/bar/bar
    aux_files.{n}/baz/baz
    """
    oldcwd = os.getcwd()
    try:
        with TemporaryDirectory() as d:
            shutil.copytree("local-props", os.path.join(d, "local-props"))
            os.chdir(d)
            os.mkdir("output")
            dotfile_path = "output/.dotfile"
            pipelinefile_path = "output/pipeline.foo"

            with open(dotfile_path, "w") as f:
                f.write("foo")
            with open(pipelinefile_path, "w") as f:
                f.write("foo")

            # Nothing happens when instantiating, only when calling
            td = FileWritingTestDriver(LennardJones())

            # Simply instantiating a second TD shouldn't break anything
            td2 = FileWritingTestDriver(LennardJones())

            # Call the TD once. Now the output directory should be established
            td()
            # Dotfile should not get touched or flagged for a nonnempty directory,
            # so nothing should have been backed up
            assert not os.path.isdir("output.0")

            assert os.path.isfile(dotfile_path)
            assert os.path.isfile(pipelinefile_path)
            assert os.path.isfile(TOKENPATH)

            # Should not go in output directory
            td.write_property_instances_to_file("results.edn")
            assert os.path.isfile("results.edn")
            n = 0
            assert os.path.isfile(f"output/foo_prop-{3*n+1}.txt")
            assert os.path.isfile(f"output/bar/bar_prop-{3*n+2}.txt")
            assert os.path.isfile(f"output/baz_prop-{3*n+3}")
            with tarfile.open(f"output/aux_files.{n}.txz") as tar:
                assert len(tar.getmembers()) == 3
                for member in tar.getmembers():
                    assert member.name in [
                        f"aux_files.{n}/foo",
                        f"aux_files.{n}/bar/bar",
                        f"aux_files.{n}/baz/baz",
                    ]
            # Baz should have been cleaned up
            assert not os.path.isdir("output/baz")

            # For checking later that after switching to a different instance,
            # output looks the same
            num_files_after_one_run = len(glob.glob("output/**"))

            # Remake baz and put a dotfile in there. Should not create a problem
            os.mkdir("output/baz")
            with open("output/baz/.dotfile", "w") as f:
                f.write("foo")

            # Run the TD again, without aux files this time
            td(write_aux_files=False)
            td.write_property_instances_to_file()
            assert os.path.isfile("output/results.edn")
            # Dotfiles should not have been touched
            assert os.path.isfile(dotfile_path)
            assert os.path.isfile(pipelinefile_path)
            assert os.path.isfile("output/baz/.dotfile")
            for n in range(2):
                assert os.path.isfile(f"output/foo_prop-{3*n+1}.txt")
                assert os.path.isfile(f"output/bar/bar_prop-{3*n+2}.txt")
                assert os.path.isfile(f"output/baz_prop-{3*n+3}")
                if n == 0:
                    with tarfile.open(f"output/aux_files.{n}.txz") as tar:
                        assert len(tar.getmembers()) == 3
                        for member in tar.getmembers():
                            assert member.name in [
                                f"aux_files.{n}/foo",
                                f"aux_files.{n}/bar/bar",
                                f"aux_files.{n}/baz/baz",
                            ]
                elif n == 1:
                    assert not os.path.exists(f"output/aux_files.{n}.txz")

            # make a non-dot file under baz. Now TD should detect the issue and crash
            with open("output/baz/nondotfile", "w") as f:
                f.write("foo")

            try:
                td()
                assert False
            except KIMTestDriverError:
                assert True

            os.mkdir("output.0")

            # Call the second Test Driver for the first time. It should find the lowest
            # output.{n} available and back up the existing output there
            td2()
            # Root dotfiles should not have been touched
            assert os.path.isfile(dotfile_path)
            assert os.path.isfile(pipelinefile_path)
            assert not os.path.exists("output.1/.dotfile")
            assert not os.path.exists("output.1/pipeline.foo")

            # Non-root dotfiles should have been moved
            assert os.path.isfile("output.1/baz/.dotfile")

            # Other files
            assert os.path.isfile("output.1/baz/nondotfile")
            assert os.path.isfile("output.1/results.edn")
            assert os.path.isfile(TOKENPATH.replace("output", "output.1"))

            # Files written inside calculate
            for n in range(2):
                assert os.path.isfile(f"output.1/foo_prop-{3*n+1}.txt")
                assert os.path.isfile(f"output.1/bar/bar_prop-{3*n+2}.txt")
                assert os.path.isfile(f"output.1/baz_prop-{3*n+3}")
                if n == 0:
                    with tarfile.open(f"output.1/aux_files.{n}.txz") as tar:
                        assert len(tar.getmembers()) == 3
                        for member in tar.getmembers():
                            assert member.name in [
                                f"aux_files.{n}/foo",
                                f"aux_files.{n}/bar/bar",
                                f"aux_files.{n}/baz/baz",
                            ]
                elif n == 1:
                    assert not os.path.exists(f"output.1/aux_files.{n}.txz")

            assert num_files_after_one_run == len(glob.glob("output/**"))

            # Try calling the first Test Driver again. Should fail
            # due to a token mismatch.
            try:
                td()
                assert False
            except KIMTestDriverError:
                assert True

    finally:
        os.chdir(oldcwd)


def test_atom_style():
    td = TestInitKIMTestDriver("LJ_ElliottAkerson_2015_Universal__MO_959249795837_003")
    assert td._get_supported_lammps_atom_style() == "atomic"
    td = TestInitKIMTestDriver(
        "Sim_LAMMPS_ReaxFF_AnGoddard_2015_BC__SM_389039364091_000"
    )
    assert td._get_supported_lammps_atom_style() == "charge"


def test_model():
    td = TestInitKIMTestDriver("LJ_ElliottAkerson_2015_Universal__MO_959249795837_003")
    assert isinstance(td.model, str)
    td = TestInitKIMTestDriver(LennardJones())
    assert isinstance(td.model, LennardJones)


if __name__ == "__main__":
    test_atom_style()
