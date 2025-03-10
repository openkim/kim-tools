#!/usr/bin/python

from kim_tools.test_driver import KIMTestDriver
from ase.atoms import Atoms
from ase.calculators.lj import LennardJones
import os

class TestTestDriver(KIMTestDriver):    
    def _calculate(self,property_name,species):
        """
        example calculate method
        
        Args:
            property_name: for testing ability to find properties at different paths.
            !!! AN ACTUAL TEST DRIVER SHOULD NOT HAVE AN ARGUMENT SUCH AS THIS !!!
        """
        atoms = Atoms([species],[[0,0,0]])
        self._add_property_instance(property_name,"This is an example disclaimer.")
        self._add_key_to_current_property_instance("species", atoms.get_chemical_symbols()[0])
        self._add_key_to_current_property_instance("mass", atoms.get_masses()[0], "amu", {'source-std-uncert-value':1})

def test_kimtest(monkeypatch):
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
        
    monkeypatch.setenv("KIM_PROPERTY_PATH", os.path.join(os.getcwd(),'mock-test-drivers-dir/*/local-props')+':'+os.path.join(os.getcwd(),'mock-test-drivers-dir/*/local_props'))

    for prop_name in testing_property_names:
        test(property_name=prop_name,species='Ar')
        
    assert len(test.property_instances) == 6
    test.write_property_instances_to_file()
    