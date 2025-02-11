#!/usr/bin/python

from kim_tools.test_driver import CrystalGenomeTestDriver, query_crystal_genome_structures
from kim_tools.aflow_util import get_stoich_reduced_list_from_prototype

class TestTestDriver(CrystalGenomeTestDriver):
    def _calculate(self,**kwargs):
        """
        example calculate method. Just writes the binding-energy and crystal-structure-npt properties assuming the provided
        structure is already at equilibrium
        """
        # calculate potential energy and do the required stuff to figure out per-formula and per-atom, and subtract isolated energy
        potential_energy = self.atoms.get_potential_energy()
        potential_energy_per_atom = potential_energy/len(self.atoms)
        reduced_stoichiometry = get_stoich_reduced_list_from_prototype(self.prototype_label) # i.e. "AB3\_...." -> [1,3]        
        binding_energy_per_formula = potential_energy_per_atom * sum(reduced_stoichiometry)
        binding_energy_per_atom = binding_energy_per_formula/sum(reduced_stoichiometry)

        # add property instance and common fields
        self._add_property_instance_and_common_crystal_genome_keys("binding-energy-crystal",False,False,"This is an example disclaimer.")

        # add the fields unique to this property
        self._add_key_to_current_property_instance("binding-potential-energy-per-atom",binding_energy_per_atom,"eV")
        self._add_key_to_current_property_instance("binding-potential-energy-per-formula",binding_energy_per_formula,"eV")

        self._add_property_instance_and_common_crystal_genome_keys("crystal-structure-npt",write_temp=True,write_stress=True)
        with open("restart.dump","w") as f:
            f.write("dummy\nfile")
        self._add_file_to_current_property_instance("restart-file","restart.dump")


def test_cg():        
    test = TestTestDriver("LJ_ElliottAkerson_2015_Universal__MO_959249795837_003")
    list_of_crystal_descriptions = query_crystal_genome_structures("LJ_ElliottAkerson_2015_Universal__MO_959249795837_003",["Ar"],"A_hP2_194_c")
    test(**list_of_crystal_descriptions[0])    
    assert len(test.property_instances) == 2
    test.write_property_instances_to_file()
    
if __name__ == '__main__':
    test_cg()