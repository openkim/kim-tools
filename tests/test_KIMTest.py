#!/usr/bin/python

from kim_tools.test_driver import KIMTestDriver
from ase.atoms import Atoms
from ase.calculators.lj import LennardJones

class TestTestDriver(KIMTestDriver):
    def _calculate(self,property_name):
        """
        example calculate method
        
        Args:
            property_name: for testing ability to find properties at different paths.
            !!! AN ACTUAL TEST DRIVER SHOULD NOT HAVE AN ARGUMENT SUCH AS THIS !!!
        """
        self.species = self.atoms.get_chemical_symbols()[0]
        self._add_property_instance(property_name,"This is an example disclaimer.")
        self._add_key_to_current_property_instance("species", self.atoms.get_chemical_symbols()[0])
        self._add_key_to_current_property_instance("mass", self.atoms.get_masses()[0], "amu", {'source-std-uncert-value':1})

def test_kimtest():
    atoms = Atoms(['Ar'], [[0, 0, 0]], cell=[[1, 0, 0], [0, 2, 0], [0, 0, 2]])
    test = TestTestDriver(LennardJones())
    testing_property_names = [
        'atomic-mass', # already in kim-properties
        'atomic-mass0', # found in $PWD/local-props
        'atomic-mass0', # check that repeat works fine
        'tag:brunnels@noreply.openkim.org,2016-05-11:property/atomic-mass1', # check that full id works as well, found in $PWD/local-props
        'atomic-mass2', # found in $PWD/local-props/atomic-mass2
        'atomic-mass3', # found in $PWD/mock-test-drivers-dir/mock-td/local_props. For testing how KDP will set these.
                        # test this using `export KIM_PROPERTY_PATH=$PWD/mock-test-drivers-dir/*/local-props:$PWD/mock-test-drivers-dir/*/local_props`
    ]

    for prop_name in testing_property_names:
        test(atoms,property_name=prop_name)
        
    assert len(test.property_instances) == 6
    test.write_property_instances_to_file()


if __name__ == '__main__':
    test_kimtest()