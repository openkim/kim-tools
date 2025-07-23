import logging
import math
import os
import shutil
from typing import Dict, List, Tuple, Union

import numpy as np
import numpy.typing as npt

from .core import DATA_DIR, _check_space_group

logger = logging.getLogger(__name__)
logging.basicConfig(filename="kim-tools.log", level=logging.INFO, force=True)


def voigt_elast_class(sgnum: Union[int, str]) -> str:
    """
    Get the name of the class of the structure of the elasticity tensor, out of
    the following possibilities:
    "cubic", "hexagonal", "trigonal_3bar_m_2nd_pos", "trigonal_3bar_m_3rd_pos",
    "trigonal_3bar", "tetragonal_class_4_slash_mmm", "tetragonal_class_4_slash_m",
    "orthorhombic", "monoclinic", "triclinic". The names that don't
    correspond to a crystal system are taken from
    https://dictionary.iucr.org/Laue_class

    Note that this is a novel classification as far as the authors are aware.
    Most crystallography texts on the subject will show 9 classes. This function
    returns 10 possibilities, because the structure of the components of
    the elasticity tensor for Laue class -3m depends on whether the twofold
    axis is aligned along the Cartesian x or y axis. Assuming one adopts a
    standard orientation of the hexagonal unit cell w.r.t. the Cartesian coordinate
    system, as we do, different space groups in this Laue class have different
    forms.

    See https://doi.org/10.1017/CBO9781139017657.008 pg 232 for a review
    of the subject.
    """
    _check_space_group(sgnum)
    sgnum = int(sgnum)

    if sgnum < 3:
        return "triclinic"
    elif sgnum < 16:
        return "monoclinic"
    elif sgnum < 75:
        return "orthorhombic"
    elif sgnum < 89:
        return "tetragonal_4_slash_m"
    elif sgnum < 143:
        return "tetragonal_4_slash_mmm"
    elif sgnum < 149:
        return "trigonal_3bar"
    elif sgnum < 168:
        # Determine if this is one of the groups with the 2-fold operation in the
        # third position (e.g. 149:P312), which has different equations
        if sgnum in (149, 151, 153, 157, 159, 162, 163):
            return "trigonal_3bar_m_3rd_pos"
        else:
            return "trigonal_3bar_m_2nd_pos"
    elif sgnum < 195:
        return "hexagonal"
    else:
        return "cubic"


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
    ELASTICITY_MATRIX_EQNS = {
        "cubic": {
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
        },
        "hexagonal": {
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
        },
        "trigonal_3bar_m_2nd_pos": {
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
        },
        "trigonal_3bar_m_3rd_pos": {
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
        },
        "trigonal_3bar": {
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
        },
        "tetragonal_4_slash_mmm": {
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
        },
        "tetragonal_4_slash_m": {
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
        },
        "orthorhombic": {
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
        },
        "monoclinic": {
            (1, 4): None,
            (1, 6): None,
            (2, 4): None,
            (2, 6): None,
            (3, 4): None,
            (3, 6): None,
            (4, 5): None,
            (5, 6): None,
        },
        "triclinic": {},
    }

    # error check typing in the above dicts
    for eqn in ELASTICITY_MATRIX_EQNS.values():
        # only unique keys
        assert sorted(list(set(eqn.keys()))) == sorted(list(eqn.keys()))
        # check that all components appearing in RHS of relations are independent, i.e.
        # they don't appear as a key
        for dependent_component in eqn:
            if eqn[dependent_component] is not None:
                for independent_component in eqn[dependent_component][1]:
                    assert not (independent_component in eqn)

    return ELASTICITY_MATRIX_EQNS[voigt_elast_class(sgnum)]


def voigt_elast_struct_svg(sgnum: Union[int, str], dest_filename: str) -> None:
    """
    Write a copy of the image showing the structure of the Voigt elasticity matrix for
    the specified space group
    """
    src_filename = os.path.join(DATA_DIR, "elast_" + voigt_elast_class(sgnum) + ".svg")
    shutil.copyfile(src_filename, dest_filename)


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
