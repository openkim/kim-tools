#!/usr/bin/python

from kim_tools.test_driver import KIMTestDriver
from ase.atoms import Atoms

class TestTestDriver(KIMTestDriver):
    def _calculate(self):
        """
        example calculate method
        """
        self.species = self.atoms.get_chemical_symbols()[0]
        self._add_property_instance("atomic-mass","This is an example disclaimer.")
        self._add_key_to_current_property_instance("species", self.atoms.get_chemical_symbols()[0])
        self._add_key_to_current_property_instance("mass", self.atoms.get_masses()[0], "amu", {'source-std-uncert-value':1})

atoms = Atoms(['Ar'], [[0, 0, 0]], cell=[[1, 0, 0], [0, 2, 0], [0, 0, 2]])
test = TestTestDriver("LennardJones_Ar")
test(atoms)
test.write_property_instances_to_file()
print(test.property_instances)
