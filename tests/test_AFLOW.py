#!/usr/bin/python

from kim_tools import AFLOW, split_parameter_array,CRYSTAL_GENOME_INITIAL_STRUCTURES, \
    get_wyckoff_lists_from_prototype, frac_pos_match_allow_permute_wrap, frac_pos_match_allow_wrap, get_real_to_virtual_species_map, \
    shuffle_atoms, get_pearson_symbol_from_prototype, minimize_wrapper, prototype_labels_are_equivalent
from ase.calculators.kim.kim import KIM
from random import random
import json
import os

TEST_CASES = [577,365,1734,1199,1478,166,1210,1362,920,212,646,22]

def test_prototype_labels_are_equivalent():
    assert not(prototype_labels_are_equivalent('AB_oC8_65_j_g','AB_oC8_65_i_i'))
    assert (prototype_labels_are_equivalent('AB2C_oP8_17_a_bc_d','AB2C_oP8_17_a_bd_c'))
    assert (prototype_labels_are_equivalent('AB_mP8_14_ad_e','AB_mP8_14_ab_e'))

def test_get_wyckoff_lists_from_prototype():
    assert get_wyckoff_lists_from_prototype('A_hP68_194_ef2h2kl') == ['efhhkkl']
    assert get_wyckoff_lists_from_prototype('AB_mC48_8_12a_12a') == ['aaaaaaaaaaaa','aaaaaaaaaaaa']
    
def test_get_param_names_from_prototype():
    aflow = AFLOW()
    assert(aflow.get_param_names_from_prototype('A_cF4_225_a') == ['a'])
    assert(aflow.get_param_names_from_prototype('A_cF240_202_h2i') == ['a','y1','z1','x2','y2','z2','x3','y3','z3'])

def test_get_equations_from_prototype():
    aflow = AFLOW()
    equation_sets_cache = {} # for large-scale testing, helpful to check that same prototype with different parameters gives the same results
    for material in [CRYSTAL_GENOME_INITIAL_STRUCTURES[test_case] for test_case in TEST_CASES]:
        species = material["species"]
        real_to_virtual_species_map = get_real_to_virtual_species_map(species)
        prototype_label = material["prototype_label"]
        parameter_names = material["parameter_names"]
        _, internal_parameter_names_ref = split_parameter_array(parameter_names)
        # TODO: Fix this
        if prototype_label.split('_')[1][:2] == 'hR':
            continue
        for parameter_set in material["parameter_sets"]:
            parameter_values = parameter_set["parameter_values"]            
            atoms = aflow.build_atoms_from_prototype(species,prototype_label,parameter_values)
            if prototype_label not in equation_sets_cache:
                equation_sets = aflow.get_equation_sets_from_prototype(prototype_label)
                equation_sets_cache[prototype_label] = equation_sets
            else:
                equation_sets = equation_sets_cache[prototype_label]
            
            internal_parameter_names = []
            for eqset in equation_sets:
                internal_parameter_names += eqset.param_names
            
            # TODO: Does it matter that the order is not the same?
            assert set(internal_parameter_names) == set(internal_parameter_names_ref), \
                f"get_equation_sets_from_prototype got incorrect internal parameter names, got {internal_parameter_names} expected {internal_parameter_names_ref}"
                
            assert sum([len(eqset.coeff_matrix_list) for eqset in equation_sets]) == len(atoms), "get_equation_sets_from_prototype got an incorrect number of equations"
            
            scaled_positions_computed_from_equations = []
            virtual_species_from_atoms = []
            species_from_equations = []
                    
            diagnostics = f'Problem occurred in prototype {prototype_label}\nReference fractional positions and species:\n'
            for atom in atoms:
                for position in atom.scaled_position:
                    diagnostics += f'{position:8.4f}'
                virtual_species = real_to_virtual_species_map[atom.symbol]
                diagnostics += f'    {virtual_species}\n'
                virtual_species_from_atoms.append(virtual_species)            
            diagnostics += '\nComputed fractional positions and species:\n'
                        
            for eqset in equation_sets:
                for coeff_mat,const_terms in zip(eqset.coeff_matrix_list,eqset.const_terms_list):
                    scaled_position = coeff_mat @ [[parameter_values[parameter_names.index(parname)]] for parname in eqset.param_names] + \
                        const_terms
                    scaled_positions_computed_from_equations.append(scaled_position.T)
                    species_from_equations.append(eqset.species)
                    for position in scaled_position:
                        diagnostics += f'{position[0]:8.4f}'
                    diagnostics += f'    {eqset.species}\n'
            
            assert frac_pos_match_allow_permute_wrap(
                atoms.get_scaled_positions(),
                scaled_positions_computed_from_equations,
                virtual_species_from_atoms,
                species_from_equations
            ), f'Failed to match fractional coordinates.\n{diagnostics}'
            
            assert frac_pos_match_allow_wrap(
                atoms.get_scaled_positions(),
                scaled_positions_computed_from_equations,
                virtual_species_from_atoms,
                species_from_equations
            ), f'Matched fractional coordinates, but there was a permutation.\n{diagnostics}'        
                        
            print(f'Successfully checked get_equations_from_prototype for label {prototype_label}')

def test_get_prototype(
    materials=[CRYSTAL_GENOME_INITIAL_STRUCTURES[test_case] for test_case in TEST_CASES]
    #materials=CRYSTAL_GENOME_INITIAL_STRUCTURES
):
    aflow = AFLOW(np=19)
    match_counts_by_pearson = {}
    match_counts_by_spacegroup = {}
    INIT_COUNTS = {'match':0,'nonmatch':0}
    
    for pearson in ['aP','mP','mC','oP','oC','oF','oI','tP','tI','hP','hR','cP','cF','cI']:
        match_counts_by_pearson[pearson] = INIT_COUNTS.copy()
    for spacegroup in range(1,231):
        match_counts_by_spacegroup[spacegroup] = INIT_COUNTS.copy()
    
    failed_to_solve_at_least_one = False
    
    for material in materials:
        species = material["species"]            
        prototype_label = material["prototype_label"]
        # TODO: Fix this
        if get_pearson_symbol_from_prototype(prototype_label).startswith('hR'):
            continue
        
        prototype_label_split = prototype_label.split('_')
        pearson = prototype_label_split[1][:2]
        spacegroup = int(prototype_label_split[2])
        
        for parameter_set in material["parameter_sets"]:
            parameter_values = parameter_set["parameter_values"]
            atoms = aflow.build_atoms_from_prototype(species,prototype_label,parameter_values)
            """
            # Using AFLOW software -- failure case for comparison
            
            cg_des = get_crystal_genome_designation_from_atoms(atoms,get_library_prototype=False,prim=False,aflow_np=19)
            
                    cg_des['stoichiometric_species'],
                    cg_des['prototype_label'],
                    cg_des['parameter_values_angstrom']
                        
            """
                
            atoms = shuffle_atoms(atoms)
            
            atoms.rotate((random(),random(),random()),(random(),random(),random()),rotate_cell=True)
            
            atoms.translate((random()*1000,random()*1000,random()*1000))
            
            atoms.wrap()
            
            atoms.calc = KIM('LJ_ElliottAkerson_2015_Universal__MO_959249795837_003')
            minimize_wrapper(atoms,fix_symmetry=True,variable_cell=True,steps=2,opt_kwargs={"maxstep":0.05})
            
            try:
                aflow.solve_for_params_of_known_prototype(atoms,prototype_label)
                crystal_did_not_rotate = True
            except Exception as e:
                crystal_did_not_rotate = False
                print(e)
            
            if not crystal_did_not_rotate:
                failed_to_solve_at_least_one = True
                print(f'Failed to solve for parameters or confirm unrotated prototype designation for {prototype_label}')
                match_counts_by_pearson[pearson]['nonmatch'] += 1
                match_counts_by_spacegroup[spacegroup]['nonmatch'] += 1
                filename = f'output/{prototype_label}.POSCAR'
                filename_exists = os.path.isfile(filename)
                if filename_exists:
                    suffix = 0
                    while not filename_exists:
                        filename = f'output/{prototype_label}.{suffix}.POSCAR'
                        filename_exists = os.path.isfile(filename)
                        suffix += 1
                print(f'Dumping atoms to {filename}')
                atoms.write(filename,format='vasp',sort=True)
            else:
                print(f'Successfully confirmed unrotated prototype designation for {prototype_label}')
                match_counts_by_pearson[pearson]['match'] += 1
                match_counts_by_spacegroup[spacegroup]['match'] += 1
        with open('output/basic_match_counts_by_pearson.json','w') as f:
            json.dump(match_counts_by_pearson,f)
        with open('output/basic_match_counts_by_spacegroup.json','w') as f:
            json.dump(match_counts_by_spacegroup,f)
            
        assert not failed_to_solve_at_least_one

if __name__ == '__main__':
    test_get_prototype()