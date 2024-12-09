#!/usr/bin/python

from kim_tools import AFLOW, split_parameter_array
import numpy as np

TEST_CASES = [
{"species": ["Ca", "O", "Si"], "prototype_label": "AB3C_aP30_2_3i_9i_3i", "parameter_names": ["a", "b/a", "c/a", "alpha", "beta", "gamma", "x1", "y1", "z1", "x2", "y2", "z2", "x3", "y3", "z3", "x4", "y4", "z4", "x5", "y5", "z5", "x6", "y6", "z6", "x7", "y7", "z7", "x8", "y8", "z8", "x9", "y9", "z9", "x10", "y10", "z10", "x11", "y11", "z11", "x12", "y12", "z12", "x13", "y13", "z13", "x14", "y14", "z14", "x15", "y15", "z15"], "parameter_values": [6.7775, 1.3825895, 0.99430468, 83.5094, 75.9601, 69.7373, 0.75709148, 0.0011147306, 0.73598857, 0.92299988, 0.33483882, 0.57277016, 0.26673272, 0.35208259, 0.90731832, 0.23178876, 0.40353274, 0.54135229, 0.61239375, 0.26425051, 0.63625742, 0.49465368, 0.13777657, 0.36197599, 0.6548869, 0.90512925, 0.11523013, 0.89834602, 0.94027959, 0.3712392, 0.80847639, 0.1347392, 0.045154102, 0.88477946, 0.40443359, 0.91223543, 0.027264404, 0.25524213, 0.23252552, 0.59809072, 0.3866145, 0.25192747, 0.47912857, 0.30143537, 0.45477333, 0.72490732, 0.015737688, 0.22911936, 0.83857573, 0.29268957, 0.10637066]},
{"species": ["C"], "prototype_label": "A_mC16_12_4i", "parameter_names": ["a", "b/a", "c/a", "beta", "x1", "z1", "x2", "z2", "x3", "z3", "x4", "z4"], "parameter_values": [9.1921, 0.27467064, 0.45133321, 82.971, 0.94271241, 0.87955035, 0.44178398, 0.65364615, 0.78567809, 0.059238811, 0.27133365, 0.58534651]},
{"species": ["U"], "prototype_label": "A_oC4_63_c", "parameter_names": ["a", "b/a", "c/a", "y1"], "parameter_values": [3.38, 1.7635799, 1.6884024, 0.14768057]},
{"species": ["Mo", "Ni"], "prototype_label": "AB4_tI10_87_a_h", "parameter_names": ["a", "c/a", "x2", "y2"], "parameter_values": [5.6931, 0.62094465, 0.40159241, 0.80168077]},
{"species": ["O", "Ti"], "prototype_label": "A2B_tP6_136_f_a", "parameter_names": ["a", "c/a", "x2"], "parameter_values": [4.6726, 0.64700595, 0.30516921]},
#{"species": ["Al", "Li"], "prototype_label": "A2B3_hR5_166_c_ac", "parameter_names": ["a", "c/a", "x2", "x3"], "parameter_values": [4.4631793, 3.1571619, 0.19703494, 0.5976076]},
#{"species": ["Al", "C"], "prototype_label": "A4B3_hR7_166_2c_ac", "parameter_names": ["a", "c/a", "x2", "x3", "x4"], "parameter_values": [3.3544538, 7.492051, 0.70650219, 0.87016163, 0.78319433]},
#{"species": ["Al", "Au"], "prototype_label": "A3B8_hR44_167_bce_2c2f", "parameter_names": ["a", "c/a", "x2", "x3", "x4", "x5", "x6", "y6", "z6", "x7", "y7", "z7"], "parameter_values": [7.8105146, 5.4347116, 0.15613556, 0.21753425, 0.93687697, 0.56692249, 0.64901599, 0.0043751172, 0.29495592, 0.61541438, 0.35227113, 0.87563146]},
#{"species": ["Al", "Pd"], "prototype_label": "AB_hR26_148_a2f_b2f", "parameter_names": ["a", "c/a", "x3", "y3", "z3", "x4", "y4", "z4", "x5", "y5", "z5", "x6", "y6", "z6"], "parameter_values": [15.732393, 0.33556832, 0.5574733, 0.84810628, 0.59879838, 0.25162269, 0.65374741, 0.096413033, 0.056906892, 0.34500705, 0.098395144, 0.75098957, 0.15501208, 0.59773116]},
{"species": ["Mo", "S"], "prototype_label": "AB2_hR3_166_a_c", "parameter_names": ["a", "c/a", "x2"], "parameter_values": [3.2047217, 5.9412636, 0.24969054]},
{"species": ["O", "Si"], "prototype_label": "A2B_hP9_154_c_a", "parameter_names": ["a", "c/a", "x1", "x2", "y2", "z2"], "parameter_values": [5.0682, 1.0946293, 0.48489843, 0.41536476, 0.23922009, 0.80786585]},
{"species": ["Fe", "P"], "prototype_label": "A2B_hP9_189_fg_ad", "parameter_names": ["a", "c/a", "x3", "x4"], "parameter_values": [6.4943, 0.52673883, 0.41418679, 0.73946286]},
{"species": ["As", "Ga"], "prototype_label": "AB_cP16_205_c_c", "parameter_names": ["a", "x1", "x2"], "parameter_values": [7.0281, 0.14318501, 0.34443472]},
{"species": ["Ce", "O"], "prototype_label": "AB2_cF12_225_a_c", "parameter_names": ["a"], "parameter_values": [5.4908]},
{"species": ["Al"], "prototype_label": "A_cF4_225_a", "parameter_names": ["a"], "parameter_values": [4.039]},
]

def test_get_equations_from_prototype():
    aflow = AFLOW()
    for case in TEST_CASES:
        parameter_names = case.pop("parameter_names")
        _, internal_parameter_names_ref = split_parameter_array(parameter_names)
        _, internal_parameter_values = split_parameter_array(parameter_names,case["parameter_values"])
        internal_parameter_values_and_one = np.asarray(internal_parameter_values+[1])
        atoms = aflow.build_atoms_from_prototype(**case)
        species = case.pop("species")
        equations, internal_parameter_names = aflow.get_equations_from_prototype(**case)

        assert internal_parameter_names == internal_parameter_names_ref, \
            f"get_equations_from_prototype got incorrect internal parameter names, got {internal_parameter_names} expected {internal_parameter_names_ref}"
        assert len(equations) == len(atoms), "get_equations_from_prototype got an incorrect number of equations"
        
        virtual_to_real_species_map = {}
        for i,symbol in enumerate(species):
            virtual_to_real_species_map[symbol]=chr(65+i)
        
        diagnostics = f'Problem occurred in prototype {case["prototype_label"]}\nReference fractional positions and species:\n'
        for atom in atoms:
            for position in atom.scaled_position:
                diagnostics += f'{position:8.4f}'
            diagnostics += f'    {virtual_to_real_species_map[atom.symbol]}\n'
        diagnostics += '\nComputed fractional positions and species:\n'
        for equation in equations:
            scaled_position = np.matmul(equation["equations"],internal_parameter_values_and_one)
            for position in scaled_position:
                diagnostics += f'{position:8.4f}'
            diagnostics += f'    {equation["species"]}\n'
            
        equation_matched = [False]*len(equations)
        for atom in atoms:
            found_match = False
            for i,equation in enumerate(equations):
                position_differences = np.matmul(equation["equations"],internal_parameter_values_and_one) - atom.scaled_position
                if np.allclose(position_differences,np.rint(position_differences),atol=1e-5):
                    assert equation["species"] == virtual_to_real_species_map[atom.symbol], \
                        f"get_equations_from_prototype got incorrect species, got {equation['species']} expected {virtual_to_real_species_map[atom.symbol]}\n\n{diagnostics}"
                    assert not found_match, f"get_equations_from_prototype gave equations that matched the same atom twice!\n\n{diagnostics}"
                    assert not equation_matched[i], f"get_equations_from_prototype gave an equations that matched two different atoms!\n\n{diagnostics}"
                    found_match = True
                    equation_matched[i] = True
            assert found_match, \
                f"get_equations_from_prototype did not produce any equations that match position {atom.scaled_position}.\n\n{diagnostics}"
                    
        print(f'Successfully checked get_equations_from_prototype for label {case["prototype_label"]}')            
        
if __name__ == '__main__':
    test_get_equations_from_prototype()