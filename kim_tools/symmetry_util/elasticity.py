import logging
import math
from typing import Dict, List, Tuple, Union

import numpy as np
import numpy.typing as npt

from .core import _check_space_group

logger = logging.getLogger(__name__)
logging.basicConfig(filename="kim-tools.log", level=logging.INFO, force=True)


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


def calc_bulk(elastic_constants):
    """
    Compute the bulk modulus given the elastic constants matrix in
    Voigt ordering.

    Parameters:
        elastic_constants : float
            A 6x6 numpy array containing the elastic constants in
            Voigt ordering. The material can have arbitrary anisotropy.

    Returns:
        bulk : float
            The bulk modulus, defined as the ratio between the hydrostatic
            stress (negative of the pressure p) in hydrostatic loading and
            the diltation e (trace of the strain tensor), i.e. B = -p/e
    """
    # Compute bulk modulus, based on exercise 6.14 in Tadmor, Miller, Elliott,
    # Continuum Mechanics and Thermodynamics, Cambridge University Press, 2012.
    rank_elastic_constants = np.linalg.matrix_rank(elastic_constants)
    elastic_constants_aug = np.concatenate(
        (elastic_constants, np.transpose([[1, 1, 1, 0, 0, 0]])), 1
    )
    rank_elastic_constants_aug = np.linalg.matrix_rank(elastic_constants_aug)
    if rank_elastic_constants_aug > rank_elastic_constants:
        assert rank_elastic_constants_aug == rank_elastic_constants + 1
        logger.info(
            "Information: Hydrostatic pressure not in the image of the elasticity "
            "matrix, zero bulk modulus!"
        )
        return 0.0
    else:
        # if a solution exists for a stress state of [1,1,1,0,0,0],
        # you can always use the pseudoinverse
        compliance = np.linalg.pinv(elastic_constants)
        bulk = 1 / np.sum(compliance[0:3, 0:3])
    return bulk


def map_to_Kelvin(C: npt.ArrayLike) -> npt.ArrayLike:
    """
    Compute the Kelvin form of the input 6x6 Voigt matrix
    """
    Ch = C.copy()
    Ch[0:3, 3:6] *= math.sqrt(2.0)
    Ch[3:6, 0:3] *= math.sqrt(2.0)
    Ch[3:6, 3:6] *= 2.0
    return Ch


def function_of_matrix(A, f):
    """Compute the function of a matrix"""
    ev, R = np.linalg.eigh(A)
    Dtilde = np.diag([f(e) for e in ev])
    return np.matmul(np.matmul(R, Dtilde), np.transpose(R))


def find_nearest_isotropy(elastic_constants):
    """
    Compute the distance between the provided matrix of elastic constants
    in Voigt notation, to the nearest matrix of elastic constants for an
    isotropic material. Return this distance, and the isotropic bulk and
    shear modulus.

    Ref: Morin, L; Gilormini, P and Derrien, K,
         "Generalized Euclidean Distances for Elasticity Tensors",
         Journal of Elasticity, Vol 138, pp. 221-232 (2020).

    Parameters:
        elastic_constants : float
            A 6x6 numpy array containing the elastic constants in
            Voigt ordering. The material can have arbitrary anisotropy.

    Returns:
        d : float
            Distance to the nearest elastic constants.
            log Euclidean metric.
        kappa : float
            Isotropic bulk modulus
        mu : float
            Isotropic shear modulus
    """
    E0 = 1.0  # arbitrary scaling constant (result unaffected by it)

    JJ = np.zeros(shape=(6, 6))
    KK = np.zeros(shape=(6, 6))
    v = {0: [0, 0], 1: [1, 1], 2: [2, 2], 3: [1, 2], 4: [0, 2], 5: [0, 1]}
    for ii in range(6):
        for jj in range(6):
            # i j k l = v[ii][0] v[ii][1] v[jj][0] v[jj][1]
            JJ[ii][jj] = (1.0 / 3.0) * (v[ii][0] == v[ii][1]) * (v[jj][0] == v[jj][1])
            KK[ii][jj] = (1.0 / 2.0) * (
                (v[ii][0] == v[jj][0]) * (v[ii][1] == v[jj][1])
                + (v[ii][0] == v[jj][1]) * (v[ii][1] == v[jj][0])
            ) - JJ[ii][jj]
    Chat = map_to_Kelvin(elastic_constants)
    JJhat = map_to_Kelvin(JJ)
    KKhat = map_to_Kelvin(KK)

    # Eqn (49) in Morin et al.
    fCoverE0 = function_of_matrix(Chat / E0, math.log)
    kappa = (E0 / 3.0) * math.exp(np.einsum("ij,ij", fCoverE0, JJhat))
    mu = (E0 / 2.0) * math.exp(0.2 * np.einsum("ij,ij", fCoverE0, KKhat))

    # Eqn (47) in Morin et al.
    dmat = (
        fCoverE0 - math.log(3.0 * kappa / E0) * JJhat - math.log(2.0 * mu / E0) * KKhat
    )
    d = math.sqrt(np.einsum("ij,ij", dmat, dmat))

    # Return results
    return d, kappa, mu
