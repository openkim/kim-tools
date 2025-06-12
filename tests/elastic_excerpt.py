# This code was originally written to symmetrize elastic
# constants matrices based on manually programming the
# algebraic relationships between components. Now we can
# use it to test the Tensor-based symmetrization

from typing import List, Tuple

import numpy as np
import numpy.typing as npt

CUBIC_EQN = {
    (1, 3): ([1], [(1, 2)]),
    (1, 4): None,
    (1, 5): None,
    (1, 6): None,
    (2, 2): ([1], [(1, 1)]),
    (2, 3): ([1], [(1, 2)]),
    (2, 4): None,
    (2, 5): None,
    (2, 6): None,
    (3, 3): ([1], [(1, 1)]),
    (3, 4): None,
    (3, 5): None,
    (3, 6): None,
    (4, 5): None,
    (4, 6): None,
    (5, 5): ([1], [(4, 4)]),
    (5, 6): None,
    (6, 6): ([1], [(4, 4)]),
}

HEXAGONAL_EQN = {
    (1, 4): None,
    (1, 5): None,
    (1, 6): None,
    (2, 2): ([1], [(1, 1)]),
    (2, 3): ([1], [(1, 3)]),
    (2, 4): None,
    (2, 5): None,
    (2, 6): None,
    (3, 4): None,
    (3, 5): None,
    (3, 6): None,
    (4, 5): None,
    (4, 6): None,
    (5, 5): ([1], [(4, 4)]),
    (5, 6): None,
    (6, 6): ([0.5, -0.5], [(1, 1), (1, 2)]),
}

TRIGONAL_CLASS_3BAR_M_SECOND_POS_EQN = {
    (1, 5): None,
    (1, 6): None,
    (2, 2): ([1], [(1, 1)]),
    (2, 3): ([1], [(1, 3)]),
    (2, 4): ([-1], [(1, 4)]),
    (2, 5): None,
    (2, 6): None,
    (3, 4): None,
    (3, 5): None,
    (3, 6): None,
    (4, 5): None,
    (4, 6): None,
    (5, 5): ([1], [(4, 4)]),
    (5, 6): ([1], [(1, 4)]),
    (6, 6): ([0.5, -0.5], [(1, 1), (1, 2)]),
}

TRIGONAL_CLASS_3BAR_M_THIRD_POS_EQN = {
    (1, 4): None,
    (1, 6): None,
    (2, 2): ([1], [(1, 1)]),
    (2, 3): ([1], [(1, 3)]),
    (2, 4): None,
    (2, 5): ([-1], [(1, 5)]),
    (2, 6): None,
    (3, 4): None,
    (3, 5): None,
    (3, 6): None,
    (4, 5): None,
    (4, 6): ([-1], [(1, 5)]),
    (5, 5): ([1], [(4, 4)]),
    (5, 6): None,
    (6, 6): ([0.5, -0.5], [(1, 1), (1, 2)]),
}

TRIGONAL_CLASS_3BAR_EQN = {
    (1, 6): None,
    (2, 2): ([1], [(1, 1)]),
    (2, 3): ([1], [(1, 3)]),
    (2, 4): ([-1], [(1, 4)]),
    (2, 5): ([-1], [(1, 5)]),
    (2, 6): None,
    (3, 4): None,
    (3, 5): None,
    (3, 6): None,
    (4, 5): None,
    (4, 6): ([-1], [(1, 5)]),
    (5, 5): ([1], [(4, 4)]),
    (5, 6): ([1], [(1, 4)]),
    (6, 6): ([0.5, -0.5], [(1, 1), (1, 2)]),
}

TETRAGONAL_CLASS_4_SLASH_MM_EQN = {
    (1, 4): None,
    (1, 5): None,
    (1, 6): None,
    (2, 2): ([1], [(1, 1)]),
    (2, 3): ([1], [(1, 3)]),
    (2, 4): None,
    (2, 5): None,
    (2, 6): None,
    (3, 4): None,
    (3, 5): None,
    (3, 6): None,
    (4, 5): None,
    (4, 6): None,
    (5, 5): ([1], [(4, 4)]),
    (5, 6): None,
}

TETRAGONAL_CLASS_4_SLASH_M_EQN = {
    (1, 4): None,
    (1, 5): None,
    (2, 2): ([1], [(1, 1)]),
    (2, 3): ([1], [(1, 3)]),
    (2, 4): None,
    (2, 5): None,
    (2, 6): ([-1], [(1, 6)]),
    (3, 4): None,
    (3, 5): None,
    (3, 6): None,
    (4, 5): None,
    (4, 6): None,
    (5, 5): ([1], [(4, 4)]),
    (5, 6): None,
}

ORTHORHOMBIC_EQN = {
    (1, 4): None,
    (1, 5): None,
    (1, 6): None,
    (2, 4): None,
    (2, 5): None,
    (2, 6): None,
    (3, 4): None,
    (3, 5): None,
    (3, 6): None,
    (4, 5): None,
    (4, 6): None,
    (5, 6): None,
}

MONOCLINIC_EQN = {
    (1, 4): None,
    (1, 6): None,
    (2, 4): None,
    (2, 6): None,
    (3, 4): None,
    (3, 6): None,
    (4, 5): None,
    (5, 6): None,
}

TRICLINIC_EQN = {}

ELASTICITY_MATRIX_EQNS = (
    CUBIC_EQN,
    HEXAGONAL_EQN,
    TRIGONAL_CLASS_3BAR_M_SECOND_POS_EQN,
    TRIGONAL_CLASS_3BAR_EQN,
    TETRAGONAL_CLASS_4_SLASH_MM_EQN,
    TETRAGONAL_CLASS_4_SLASH_M_EQN,
    ORTHORHOMBIC_EQN,
    MONOCLINIC_EQN,
    TRICLINIC_EQN,
)

# error check typing in the above dicts
for eqn in ELASTICITY_MATRIX_EQNS:
    assert sorted(list(set(eqn.keys()))) == sorted(list(eqn.keys()))  # only unique keys
    # check that all components appearing in RHS of relations are independent, i.e. they
    # don't appear as a key
    for dependent_component in eqn:
        if eqn[dependent_component] is not None:
            for independent_component in eqn[dependent_component][1]:
                assert not (independent_component in eqn)

# TODO: Rework this


def get_unique_components_and_reconstruct_matrix(
    elastic_constants: npt.ArrayLike, space_group_number: int
) -> Tuple[List[str], List[float], npt.ArrayLike]:
    """
    From an elasticity matrix in Voigt order and a space group number, extract the
    elastic constants that should be unique (cij where first i is as low as possible,
    then j). Reconstruct the elasticity matrix based on the algebraic symmetry rules
    and check how much the original matrix violated the symmetry rules for both
    crystallography and material frame indifference

    Returns:
        * Names of unique elastic constants
        * List of unique elastic constants
        * Reconstructed matrix
    """
    assert 0 < space_group_number < 231

    if space_group_number < 3:
        eqn = TRICLINIC_EQN
    elif space_group_number < 16:
        eqn = MONOCLINIC_EQN
    elif space_group_number < 75:
        eqn = ORTHORHOMBIC_EQN
    elif space_group_number < 89:
        eqn = TETRAGONAL_CLASS_4_SLASH_M_EQN
    elif space_group_number < 143:
        eqn = TETRAGONAL_CLASS_4_SLASH_MM_EQN
    elif space_group_number < 149:
        eqn = TRIGONAL_CLASS_3BAR_EQN
    elif space_group_number < 168:
        eqn = TRIGONAL_CLASS_3BAR_M_SECOND_POS_EQN
    elif space_group_number < 195:
        eqn = HEXAGONAL_EQN
    else:
        eqn = CUBIC_EQN

    if eqn == TRIGONAL_CLASS_3BAR_M_SECOND_POS_EQN:
        # Determine if this is one of the groups with the 2-fold operation in the
        # third position (e.g. 149:P312), which has different equations
        if space_group_number in (149, 151, 153, 157, 159, 162, 163):
            eqn = TRIGONAL_CLASS_3BAR_M_THIRD_POS_EQN

    elastic_constants_names = []
    elastic_constants_values = []
    reconstructed_matrix = np.zeros((6, 6))
    # first, figure out which constants are unique and extract them
    for i in range(1, 7):
        for j in range(i, 7):
            if (i, j) not in eqn:
                elastic_constants_names.append("c" + str(i) + str(j))
                elastic_constants_values.append(elastic_constants[i - 1, j - 1])
                reconstructed_matrix[i - 1, j - 1] = elastic_constants[i - 1, j - 1]

    for dep_comp in eqn:
        if eqn[dep_comp] is None:
            continue
        else:
            dep_comp_zero_based = (dep_comp[0] - 1, dep_comp[1] - 1)
            for coeff, indep_comp in zip(*eqn[dep_comp]):
                indep_comp_zero_based = (indep_comp[0] - 1, indep_comp[1] - 1)
                reconstructed_matrix[dep_comp_zero_based] += (
                    coeff * reconstructed_matrix[indep_comp_zero_based]
                )

    reconstructed_matrix = (
        reconstructed_matrix
        + reconstructed_matrix.T
        - np.diag(reconstructed_matrix.diagonal())
    )

    return elastic_constants_names, elastic_constants_values, reconstructed_matrix
