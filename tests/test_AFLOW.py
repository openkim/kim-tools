from kim_tools import AFLOW,CRYSTAL_GENOME_INITIAL_STRUCTURES
from ase.spacegroup import Spacegroup
import json
import numpy as np
from scipy.spatial.transform import Rotation 
proto_poscar = 'output/tmp_proto.POSCAR'
spglib_poscar = 'output/tmp_spglib.POSCAR'
output_file = 'output/grades.json'
aflow = AFLOW(np=19)

results = []

for structure in CRYSTAL_GENOME_INITIAL_STRUCTURES:
    species = structure['species']
    prototype_label = structure['prototype_label']
    for parameter_set in structure['parameter_sets']:            
        parameter_values = parameter_set['parameter_values']
        aflow_designation={'species':species,'prototype_label':prototype_label,'parameter_values':parameter_values}

        print()
        print('aflow --proto=%s:%s --params=%s'%(prototype_label,':'.join(species),','.join([str(value) for value in parameter_values])))
        print()
        try:
            atoms = aflow.build_atoms_from_prototype(**aflow_designation,proto_file=proto_poscar,verbose=False)
            atoms.write(spglib_poscar,format='vasp',sort=True)
            space_group = int(prototype_label.split('_')[2])

            _,cart_rot,frac_trans_with_ones = aflow.get_basistransformation_rotation_originshift_from_poscars(spglib_poscar,proto_poscar)

            print('\nTRANSFORMATION FROM --proto TO SPGLIB:')

            print('\nROTATION (CARTESIAN):\n')
            print(np.array_str(cart_rot,suppress_small=True))

            frac_trans = []
            for component in frac_trans_with_ones:
                if np.allclose(component,-1.0) or np.allclose(component,1.0): # sometimes can happen in aflow for some reason
                    frac_trans.append(0.)
                else:
                    frac_trans.append(component)

            frac_trans = np.asarray(frac_trans)

            if prototype_label.split('_')[1][0]=='h':
                # we expect to have a 60 degree rotation going from proto to spglib, so undo it
                cart_rot = np.matmul(cart_rot,Rotation.from_rotvec([0,0,-60],True).as_matrix())
                print('\nROTATION ACCOUNTING FOR HEXAGONAL CELL CHANGE (CARTESIAN):\n')
                print(np.array_str(cart_rot,suppress_small=True))
                
            # I think using the transformed cell is fine here, since we are assuming the cells are identical
            cryst_rot = np.matmul(np.linalg.inv(atoms.cell.T),np.matmul(cart_rot,atoms.cell.T))

            print('\nROTATION (CELL):\n')
            print(np.array_str(cryst_rot,suppress_small=True))

            print('\nTRANSLATION (CELL):\n')
            print(np.array_str(frac_trans,suppress_small=True))

            found_match = False
            for symop in Spacegroup(space_group).get_symop():
                if np.allclose(cryst_rot,symop[0],atol=1e-4) and np.allclose(frac_trans,symop[1],atol=1e-4):
                    print('\n!!! PASS, OPERATION FOUND IN SPACE GROUP %d !!!\n'%space_group)
                    results.append({'Unchanged':True})
                    found_match = True
                    break
            if not found_match:
                print('\n!!! FAIL, OPERATION NOT FOUND IN SPACE GROUP %d !!!\n'%space_group)
                results.append({'Unchanged':False})

            results[-1].update(aflow_designation)
        except Exception as e:
            print('\n!!! EXCEPTION ENCOUNTERED: !!!\n')
            print()
            print(e)
            print()

with open(output_file,'w') as f:
    json.dump(results,f)