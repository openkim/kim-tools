from ase.spacegroup import Spacegroup,get_spacegroup
import numpy as np

#def cartesian_rotation_is_in_space_group(rotation_matrix,cell,space_group) -> bool:
#    # if space_group is None:  # add detection      
#    for symop in Spacegroup(space_group).get_symop():
#        if

def cryst_rot_from_cart_rot(cart_rot,cell):
    return np.matmul(np.linalg.inv(cell.T),np.matmul(cart_rot,cell.T))