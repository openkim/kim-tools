#!/usr/bin/python

from kim_test_utils.test_driver import CrystalGenomeTestDriver, query_crystal_genome_structures
from kim_python_utils.ase import get_isolated_energy_per_atom
from crystal_genome_util.aflow_util import get_stoich_reduced_list_from_prototype

class TestTestDriver(CrystalGenomeTestDriver):
    def _calculate(self,**kwargs):
        """
        example calculate method. Just writes the binding-energy and crystal-structure-npt properties assuming the provided
        structure is already at equilibrium

        Args:
            structure_index:
                KIM tests can loop over multiple structures (i.e. crystals, molecules, etc.). This indicates which is being used for the current calculation.        
        """
        # calculate potential energy and do the required stuff to figure out per-formula and per-atom, and subtract isolated energy
        potential_energy = self.atoms.get_potential_energy()
        potential_energy_per_atom = potential_energy/len(self.atoms)
        reduced_stoichiometry = get_stoich_reduced_list_from_prototype(self.prototype_label) # i.e. "AB3\_...." -> [1,3]        
        binding_energy_per_formula = potential_energy_per_atom * sum(reduced_stoichiometry)
        for num_in_formula,species in zip(reduced_stoichiometry,self.stoichiometric_species):
            binding_energy_per_formula -= num_in_formula*get_isolated_energy_per_atom(self.kim_model_name,species)
        binding_energy_per_atom = binding_energy_per_formula/sum(reduced_stoichiometry)

        # add property instance and common fields
        self._add_property_instance_and_common_crystal_genome_keys("binding-energy-crystal",write_stress=False, write_temp=False)

        # add the fields unique to this property
        self._add_key_to_current_property_instance("binding-potential-energy-per-atom",binding_energy_per_atom,"eV")
        self._add_key_to_current_property_instance("binding-potential-energy-per-formula",binding_energy_per_formula,"eV")

        
test = TestTestDriver("MEAM_LAMMPS_KoJimLee_2012_FeP__MO_179420363944_002")
list_of_crystal_descriptions = query_crystal_genome_structures("MEAM_LAMMPS_KoJimLee_2012_FeP__MO_179420363944_002",["Fe","P"],"AB_oP8_62_c_c")
print(list_of_crystal_descriptions)
test(**list_of_crystal_descriptions[0])
test(**list_of_crystal_descriptions[1])
print(test.get_property_instances())
test.write_property_instances_to_file()