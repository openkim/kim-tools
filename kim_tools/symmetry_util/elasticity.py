from typing import Dict, List, Tuple, Union

import numpy.typing as npt

from .core import _check_space_group


def voigt_elast_compon_eqn(sgnum: Union[int, str]) -> Dict:
    """
    Get the algebraic equations describing the symmetry restrictions
    on the elasticity matrix in Voigt form given a space group number.
    The unit cell
    must be in the orientation defined in doi.org/10.1016/j.commatsci.2017.01.017
    for these equations to be correct.

    Returns:
        Encoding of symmetry restrictions on elasticity matrices.
        The keys are the Voigt indices of non-independent components.
        The values are a pair of lists representing the linear combination
        of the unique compoonents that is used to determine the non-unique component
        specified in the key. The first list is the coefficients, the second
        is the indices. If a non-independent component is zero, this is indicated
        by a value of None. Any components not listed as a key are assumed to
        be independent. Only the upper triangle (i<j) is listed. Indices are
        one-based.
    """
    _check_space_group(sgnum)
    sgnum = int(sgnum)

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
        TRIGONAL_CLASS_3BAR_M_THIRD_POS_EQN,
        TRIGONAL_CLASS_3BAR_EQN,
        TETRAGONAL_CLASS_4_SLASH_MM_EQN,
        TETRAGONAL_CLASS_4_SLASH_M_EQN,
        ORTHORHOMBIC_EQN,
        MONOCLINIC_EQN,
        TRICLINIC_EQN,
    )

    # error check typing in the above dicts
    for eqn in ELASTICITY_MATRIX_EQNS:
        # only unique keys
        assert sorted(list(set(eqn.keys()))) == sorted(list(eqn.keys()))
        # check that all components appearing in RHS of relations are independent, i.e.
        # they don't appear as a key
        for dependent_component in eqn:
            if eqn[dependent_component] is not None:
                for independent_component in eqn[dependent_component][1]:
                    assert not (independent_component in eqn)

    if sgnum < 3:
        eqn = TRICLINIC_EQN
    elif sgnum < 16:
        eqn = MONOCLINIC_EQN
    elif sgnum < 75:
        eqn = ORTHORHOMBIC_EQN
    elif sgnum < 89:
        eqn = TETRAGONAL_CLASS_4_SLASH_M_EQN
    elif sgnum < 143:
        eqn = TETRAGONAL_CLASS_4_SLASH_MM_EQN
    elif sgnum < 149:
        eqn = TRIGONAL_CLASS_3BAR_EQN
    elif sgnum < 168:
        eqn = TRIGONAL_CLASS_3BAR_M_SECOND_POS_EQN
    elif sgnum < 195:
        eqn = HEXAGONAL_EQN
    else:
        eqn = CUBIC_EQN

    if eqn == TRIGONAL_CLASS_3BAR_M_SECOND_POS_EQN:
        # Determine if this is one of the groups with the 2-fold operation in the
        # third position (e.g. 149:P312), which has different equations
        if sgnum in (149, 151, 153, 157, 159, 162, 163):
            eqn = TRIGONAL_CLASS_3BAR_M_THIRD_POS_EQN

    return eqn


def indep_elast_compon_names_and_values_from_voigt(
    voigt: npt.ArrayLike, sgnum: Union[int, str]
) -> Tuple[List[str], List[float]]:
    """
    From an elasticity matrix in Voigt order and a space group number,
    extract the elastic constants that should be unique (cij where first i is as low as
    possible, then j)
    """
    eqn = voigt_elast_compon_eqn(sgnum)

    elastic_constants_names = []
    elastic_constants_values = []

    # first, figure out which constants are unique and extract them
    for i in range(1, 7):
        for j in range(i, 7):
            if (i, j) not in eqn:
                elastic_constants_names.append("c" + str(i) + str(j))
                elastic_constants_values.append(voigt[i - 1, j - 1])

    return elastic_constants_names, elastic_constants_values
