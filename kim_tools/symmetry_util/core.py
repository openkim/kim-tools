"""
Crystal Symmetry utilities and data that are (mostly) independent of AFLOW
"""

import json
import logging
import os
from itertools import product
from math import ceil
from typing import Dict, List, Union

import numpy as np
from ase import Atoms
from ase.cell import Cell
from ase.geometry import get_duplicate_atoms
from numpy.typing import ArrayLike
from sympy import Matrix, cos, matrix2numpy, sin, sqrt, symbols

logger = logging.getLogger(__name__)
logging.basicConfig(filename="kim-tools.log", level=logging.INFO, force=True)

__all__ = [
    "BRAVAIS_LATTICES",
    "FORMAL_BRAVAIS_LATTICES",
    "CENTERING_DIVISORS",
    "C_CENTERED_ORTHORHOMBIC_GROUPS",
    "A_CENTERED_ORTHORHOMBIC_GROUPS",
    "IncorrectCrystallographyException",
    "IncorrectNumAtomsException",
    "are_in_same_wyckoff_set",
    "space_group_numbers_are_enantiomorphic",
    "cartesian_to_fractional_itc_rotation_from_ase_cell",
    "cartesian_rotation_is_in_point_group",
    "get_cell_from_poscar",
    "get_wyck_pos_xform_under_normalizer",
    "get_bravais_lattice_from_space_group",
    "get_formal_bravais_lattice_from_space_group",
    "get_primitive_wyckoff_multiplicity",
    "get_symbolic_cell_from_formal_bravais_lattice",
    "get_change_of_basis_matrix_to_conventional_cell_from_formal_bravais_lattice",
    "change_of_basis_atoms",
    "get_possible_primitive_shifts",
    "get_primitive_genpos_ops",
]

C_CENTERED_ORTHORHOMBIC_GROUPS = (20, 21, 35, 36, 37, 63, 64, 65, 66, 67, 68)
A_CENTERED_ORTHORHOMBIC_GROUPS = (38, 39, 40, 41)
BRAVAIS_LATTICES = [
    "aP",
    "mP",
    "mC",
    "oP",
    "oC",
    "oI",
    "oF",
    "tP",
    "tI",
    "hP",
    "hR",
    "cP",
    "cF",
    "cI",
]
FORMAL_BRAVAIS_LATTICES = BRAVAIS_LATTICES + ["oA"]
CENTERING_DIVISORS = {
    "P": 1,
    "C": 2,
    "A": 2,
    "I": 2,
    "F": 4,
    "R": 3,
}
DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")


class IncorrectCrystallographyException(Exception):
    """
    Raised when incorrect data is provided, e.g. nonexistent Bravais lattice etc.
    """


class IncorrectNumAtomsException(Exception):
    """
    Raised when the a disagreement in the number of atoms is found.
    """


def _check_space_group(sgnum: Union[int, str]):
    try:
        assert 1 <= int(sgnum) <= 230
    except Exception:
        raise IncorrectCrystallographyException(
            f"Got a space group number {sgnum} that is non-numeric or not between 1 "
            "and 230 inclusive"
        )


def cartesian_to_fractional_itc_rotation_from_ase_cell(
    cart_rot: ArrayLike, cell: ArrayLike
) -> ArrayLike:
    """
    Convert Cartesian to fractional rotation. Read the arguments and returns carefully,
    as there is some unfortunate mixing of row and columns because of the different
    conventions of the ITC and ASE and other simulation packages

    Args:
        cart_rot:
            Cartesian rotation. It is assumed that this is for left-multiplying column
            vectors, although in cases where we don't care if we're working with the
            rotation or its inverse (e.g. when checking whether or not it's in the
            point group), this doesn't matter due to orthogonality
        cell:
            The cell of the crystal, with each row being a cartesian vector
            representing a lattice vector. This is consistent with most simulation
            packages, but transposed from the ITC

    Returns:
        The fractional rotation in ITC convention, i.e. for left-multiplying column
        vectors. Here the distinction with a matrix's transpose DOES matter, because
        the fractional coordinate system is not orthonormal.
    """

    cell_arr = np.asarray(cell)
    cart_rot_arr = np.asarray(cart_rot)

    if not ((cell_arr.shape == (3, 3)) and (cart_rot_arr.shape == (3, 3))):
        raise IncorrectCrystallographyException(
            "Either the rotation matrix or the cell provided were not 3x3 matrices"
        )

    return np.transpose(cell_arr @ cart_rot_arr @ np.linalg.inv(cell_arr))


def fractional_to_cartesian_itc_rotation_from_ase_cell(
    frac_rot: ArrayLike, cell: ArrayLike
) -> ArrayLike:
    """
    Convert fractional to Cartesian rotation. Read the arguments and returns carefully,
    as there is some unfortunate mixing of row and columns because of the different
    conventions of the ITC and ASE and other simulation packages

    Args:
        frac_rot:
            The fractional rotation in ITC convention, i.e. for left-multiplying column
            vectors. Here the distinction with a matrix's transpose DOES matter, because
            the fractional coordinate system is not orthonormal.
        cell:
            The cell of the crystal, with each row being a cartesian vector
            representing a lattice vector. This is consistent with most simulation
            packages, but transposed from the ITC

    Returns:
        Cartesian rotation. It is assumed that this is for left-multiplying column
        vectors, although in cases where we don't care if we're working with the
        rotation or its inverse (e.g. when checking whether or not it's in the
        point group), this doesn't matter due to orthogonality
    """

    cell_arr = np.asarray(cell)
    frac_rot_arr = np.asarray(frac_rot)

    if not ((cell_arr.shape == (3, 3)) and (frac_rot_arr.shape == (3, 3))):
        raise IncorrectCrystallographyException(
            "Either the rotation matrix or the cell provided were not 3x3 matrices"
        )

    return np.transpose(np.linalg.inv(cell_arr) @ frac_rot_arr @ cell_arr)


def cartesian_rotation_is_in_point_group(
    cart_rot: ArrayLike, sgnum: Union[int, str], cell: ArrayLike
) -> bool:
    """
    Check that a Cartesian rotation is in the point group of a crystal given by its
    space group number and primitive cell

    Args:
        cart_rot:
            Cartesian rotation
        sgnum:
            space group number
        cell:
            The *primitive* cell of the crystal as defined in
            http://doi.org/10.1016/j.commatsci.2017.01.017, with each row being a
            cartesian vector representing a lattice vector. This is
            consistent with most simulation packages, but transposed from the ITC
    """
    # we don't care about properly transposing (i.e. worrying whether it's operating on
    # row or column vectors) the input cart_rot because that one is orthogonal, and
    # both it and its inverse must be in the point group
    frac_rot = cartesian_to_fractional_itc_rotation_from_ase_cell(cart_rot, cell)

    space_group_ops = get_primitive_genpos_ops(sgnum)

    for op in space_group_ops:
        if np.allclose(frac_rot, op["W"], atol=1e-4):
            logger.info("Found matching rotation")
            return True

    logger.info("No matching rotation found")
    return False


def get_cell_from_poscar(poscar: os.PathLike) -> ArrayLike:
    """
    Extract the unit cell from a POSCAR file, including the specified scaling
    """
    with open(poscar) as f:
        poscar_lines = f.read().splitlines()

    scaling = float(poscar_lines[1])
    cell = np.asarray(
        [[float(num) for num in line.split()] for line in poscar_lines[2:5]]
    )

    if scaling < 0:
        desired_volume = -scaling
        unscaled_volume = Cell(cell).volume
        scaling = (desired_volume / unscaled_volume) ** (1 / 3)

    return cell * scaling


def are_in_same_wyckoff_set(letter_1: str, letter_2: str, sgnum: Union[str, int]):
    """
    Given two Wyckoff letters and a space group number, return whether or not they are
    in the same Wyckoff set, meaning that their orbits are related by an operation in
    the normalizer of the space group
    """
    _check_space_group(sgnum)
    with open(os.path.join(DATA_DIR, "wyckoff_sets.json")) as f:
        wyckoff_sets = json.load(f)
    for wyckoff_set in wyckoff_sets[str(sgnum)]:
        if letter_1 in wyckoff_set:
            if letter_2 in wyckoff_set:
                return True
            else:
                return False


def space_group_numbers_are_enantiomorphic(sg_1: int, sg_2: int) -> bool:
    """
    Return whether or not two spacegroups (specified by number) are enantiomorphs of
    each other
    """
    _check_space_group(sg_1)
    _check_space_group(sg_2)
    if sg_1 == sg_2:
        return True
    else:
        enantiomorph_conversion = {
            78: 76,
            95: 91,
            96: 92,
            145: 144,
            153: 151,
            154: 152,
            170: 169,
            172: 171,
            179: 178,
            181: 180,
            213: 212,
        }
        enantiomorph_conversion_2 = {v: k for k, v in enantiomorph_conversion.items()}
        enantiomorph_conversion.update(enantiomorph_conversion_2)
        if enantiomorph_conversion[sg_1] == sg_2:
            return True
        else:
            return False


def get_wyck_pos_xform_under_normalizer(sgnum: Union[int, str]) -> List[List[str]]:
    """
    Get the "Transformed WP" column of the tables at the bottom of the page for each
    space group from https://cryst.ehu.es/cryst/get_set.html
    """
    _check_space_group(sgnum)
    with open(os.path.join(DATA_DIR, "wyck_pos_xform_under_normalizer.json")) as f:
        wyck_pos_xform_under_normalizer = json.load(f)
    return wyck_pos_xform_under_normalizer[str(sgnum)]


def get_bravais_lattice_from_space_group(sgnum: Union[int, str]):
    """
    Get the symbol (e.g. 'cF') of one of the 14 Bravais lattices from the space group
    number
    """
    _check_space_group(sgnum)
    with open(
        os.path.join(DATA_DIR, "space_groups_for_each_bravais_lattice.json")
    ) as f:
        space_groups_for_each_bravais_lattice = json.load(f)
    for bravais_lattice in space_groups_for_each_bravais_lattice:
        if int(sgnum) in space_groups_for_each_bravais_lattice[bravais_lattice]:
            return bravais_lattice
    raise RuntimeError(
        f"Failed to find space group number f{sgnum} in table of lattice symbols"
    )


def get_formal_bravais_lattice_from_space_group(sgnum: Union[int, str]):
    """
    Same as :func:`get_bravais_lattice_from_space_group` except distinguish between "oA"
    and "oC"
    """
    bravais_lattice = get_bravais_lattice_from_space_group(sgnum)
    if bravais_lattice == "oC":
        if int(sgnum) in A_CENTERED_ORTHORHOMBIC_GROUPS:
            return "oA"
        else:
            assert int(sgnum) in C_CENTERED_ORTHORHOMBIC_GROUPS
    return bravais_lattice


def get_primitive_wyckoff_multiplicity(sgnum: Union[int, str], wyckoff: str) -> int:
    """
    Get the multiplicity of a given Wyckoff letter for a primitive cell of the crystal
    """
    _check_space_group(sgnum)
    centering_divisor = CENTERING_DIVISORS[
        get_bravais_lattice_from_space_group(sgnum)[1]
    ]
    with open(os.path.join(DATA_DIR, "wyckoff_multiplicities.json")) as f:
        wyckoff_multiplicities = json.load(f)
    multiplicity_per_primitive_cell = (
        wyckoff_multiplicities[str(sgnum)][wyckoff] / centering_divisor
    )
    # check that multiplicity is an integer
    assert np.isclose(
        multiplicity_per_primitive_cell, round(multiplicity_per_primitive_cell)
    )
    return round(multiplicity_per_primitive_cell)


def get_symbolic_cell_from_formal_bravais_lattice(
    formal_bravais_lattice: str,
) -> Matrix:
    """
    Get the symbolic primitive unit cell as defined in
    http://doi.org/10.1016/j.commatsci.2017.01.017 in terms of the appropriate
    (possibly trivial) subset of the parameters a, b, c, alpha, beta, gamma

    Args:
        formal_bravais_lattice:
            The symbol for the Bravais lattice, e.g "oA". Specifically, "oA" is
            distinguished from "oC", meaning there are 15 possibilities, not just the
            14 Bravais lattices.

    Returns:
        Symbolic 3x3 matrix with the rows being cell vectors. This is in agreement with
        most simulation software, but the transpose of how the ITA defines cell vectors.

    Raises:
        IncorrectCrystallographyException:
            If a nonexistent Bravais lattice is provided
    """
    if formal_bravais_lattice not in FORMAL_BRAVAIS_LATTICES:
        raise IncorrectCrystallographyException(
            f"The provided Bravais lattice type {formal_bravais_lattice} "
            "does not exist."
        )

    a, b, c, alpha, beta, gamma = symbols("a b c alpha beta gamma")

    if formal_bravais_lattice == "aP":
        c_x = c * cos(beta)
        c_y = c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma)
        c_z = sqrt(c**2 - c_x**2 - c_y**2)
        return Matrix([[a, 0, 0], [0, b * cos(gamma), b * sin(gamma)], [c_x, c_y, c_z]])
    elif formal_bravais_lattice == "mP":
        return Matrix([[a, 0, 0], [0, b, 0], [c * cos(beta), 0, c * sin(beta)]])
    elif formal_bravais_lattice == "mC":
        return Matrix(
            [[a / 2, -b / 2, 0], [a / 2, b / 2, 0], [c * cos(beta), 0, c * sin(beta)]]
        )
    elif formal_bravais_lattice == "oP":
        return Matrix([[a, 0, 0], [0, b, 0], [0, 0, c]])
    elif formal_bravais_lattice == "oC":
        return Matrix([[a / 2, -b / 2, 0], [a / 2, b / 2, 0], [0, 0, c]])
    elif formal_bravais_lattice == "oA":
        return Matrix([[a, 0, 0], [0, b / 2, -c / 2], [0, b / 2, c / 2]])
    elif formal_bravais_lattice == "oI":
        return Matrix(
            [[-a / 2, b / 2, c / 2], [a / 2, -b / 2, c / 2], [a / 2, b / 2, -c / 2]]
        )
    elif formal_bravais_lattice == "oF":
        return Matrix([[0, b / 2, c / 2], [a / 2, 0, c / 2], [a / 2, b / 2, 0]])
    elif formal_bravais_lattice == "tP":
        return Matrix([[a, 0, 0], [0, a, 0], [0, 0, c]])
    elif formal_bravais_lattice == "tI":
        return Matrix(
            [[-a / 2, a / 2, c / 2], [a / 2, -a / 2, c / 2], [a / 2, a / 2, -c / 2]]
        )
    elif formal_bravais_lattice == "hP":
        return Matrix(
            [[a / 2, -sqrt(3) * a / 2, 0], [a / 2, sqrt(3) * a / 2, 0], [0, 0, c]]
        )
    elif formal_bravais_lattice == "hR":
        return Matrix(
            [
                [a / 2, -a / (2 * sqrt(3)), c / 3],
                [0, a / sqrt(3), c / 3],
                [-a / 2, -a / (2 * sqrt(3)), c / 3],
            ]
        )
    elif formal_bravais_lattice == "cP":
        return Matrix(
            [
                [a, 0, 0],
                [0, a, 0],
                [0, 0, a],
            ]
        )
    elif formal_bravais_lattice == "cI":
        return Matrix(
            [[-a / 2, a / 2, a / 2], [a / 2, -a / 2, a / 2], [a / 2, a / 2, -a / 2]]
        )
    elif formal_bravais_lattice == "cF":
        return Matrix([[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]])
    else:
        assert False


def get_change_of_basis_matrix_to_conventional_cell_from_formal_bravais_lattice(
    formal_bravais_lattice: str,
) -> ArrayLike:
    """
    Get a change of basis matrix **P** as defined in ITA 1.5.1.2, with "old basis"
    being the primitive cell of the provided Bravais lattice, and the "new basis" being
    the conventional cell, i.e. the cell of the primitive lattice of the same crystal
    family. E.g. if ``formal_bravais_lattice="oA"``, then "old basis" is oA, and
    "new basis" is oP. The cell choices are defined in
    http://doi.org/10.1016/j.commatsci.2017.01.017,
    including distinguishing between oA and oC.

    The matrices are given in ITC convention, meaning that they expect to operate on
    column vectors, i.e. The bases are related by the following, where the primed
    symbols indicate the new basis:

    Relationship between basis vectors:
    (**a**', **b**', **c**') = (**a**, **b**, **c**) **P**

    Relationship between fractional coordinates in each basis: **x** = **P** **x**'

    For operating on row vectors, as is often given in simulation software, make sure
    to transpose these relationships appropriately.

    Args:
        formal_bravais_lattice:
            The symbol for the Bravais lattice, e.g "oA". Specifically, "oA" is
            distinguished from "oC", meaning there are 15 possibilities, not just the
            14 Bravais lattices.

    Returns:
        Integral 3x3 matrix representing the change of basis

    Raises:
        IncorrectCrystallographyException:
            If a nonexistent Bravais lattice is provided
    """
    if formal_bravais_lattice not in FORMAL_BRAVAIS_LATTICES:
        raise IncorrectCrystallographyException(
            f"The provided Bravais lattice type {formal_bravais_lattice} "
            "does not exist."
        )

    if formal_bravais_lattice[1] == "P":  # Already primitive
        return np.eye(3)

    corresponding_primitive_lattice = formal_bravais_lattice[0] + "P"

    old_basis = get_symbolic_cell_from_formal_bravais_lattice(
        formal_bravais_lattice
    ).transpose()
    new_basis = get_symbolic_cell_from_formal_bravais_lattice(
        corresponding_primitive_lattice
    ).transpose()

    change_of_basis_matrix = matrix2numpy((old_basis**-1) @ new_basis, dtype=float)

    # matrices should be integral
    assert np.allclose(np.round(change_of_basis_matrix), change_of_basis_matrix)

    return np.round(change_of_basis_matrix)


def change_of_basis_atoms(atoms: Atoms, change_of_basis: ArrayLike) -> Atoms:
    """
    Perform an arbitrary basis change on an ``Atoms`` object, duplicating or cropping
    atoms as needed. A basic check is made that the determinant of ``change_of_basis``
    is compatible with the number of atoms, but this is far from fully determining
    that ``change_of_basis`` is appropriate for the particuar crystal described by
    ``atoms``, which is up to the user.

    NOTE: This requires ASE >= 3.25 to delete atoms that are close across PBCs

    Args:
        atoms:
            The object to transform
        change_of_basis:
            A change of basis matrix **P** as defined in ITA 1.5.1.2, with ``atoms``
            corresponding to the "old basis" and the returned ``Atoms`` object being
            in the "new basis".

            This matrix should be given in ITC convention, meaning that it expects to
            operate on column vectors, i.e. The bases are related by the following,
            where the primed symbols indicate the new basis:

            Relationship between basis vectors:
            (**a**', **b**', **c**') = (**a**, **b**, **c**) **P**

            Relationship between fractional coordinates in each basis:
            **x** = **P** **x**'

    Returns:
        The transformed ``Atoms`` object, containing the original number of
        atoms mutiplied by the determinant of the change of basis.
    """
    old_cell_column = np.transpose(atoms.cell)
    new_cell_column = old_cell_column @ change_of_basis
    new_cell = np.transpose(new_cell_column)

    # There are surely better ways to do this, but the simplest way I can think of
    # is simply to use ``Atoms.repeat()`` to create a supercell big enough to encase
    # the ``new_cell``, then wrap the atoms back into ``new_cell`` and delete dupes
    repeat = []
    for old_cell_vector in atoms.cell:
        this_repeat = 0
        old_cell_vector_norm = np.linalg.norm(old_cell_vector)
        old_cell_vector_unit = old_cell_vector / old_cell_vector_norm
        # We need to repeat the old vector enough times that it is big enough
        # to cover all possible combinations of projected new vectors
        projections = [
            np.dot(new_cell_vector, old_cell_vector_unit)
            for new_cell_vector in new_cell
        ]
        absmax_projected_sum = 0
        for coeffs in product((-1, 1), repeat=3):
            projected_sum = np.dot(coeffs, projections)
            absmax_projected_sum = max(absmax_projected_sum, abs(projected_sum))
        absmax_projected_sum += 0.1  # pad it a little bit
        this_repeat = ceil(absmax_projected_sum / old_cell_vector_norm)
        repeat.append(this_repeat)

    new_atoms = atoms.repeat(repeat)
    new_atoms.set_cell(new_cell)
    new_atoms.wrap()
    get_duplicate_atoms(new_atoms, delete=True)

    volume_change = np.linalg.det(change_of_basis)
    if not np.isclose(len(atoms) * volume_change, len(new_atoms)):
        raise IncorrectNumAtomsException(
            f"The change in the number of atoms from {len(atoms)} to {len(new_atoms)} "
            f"disagrees with the fractional change in cell volume {volume_change}"
        )

    return new_atoms


def get_possible_primitive_shifts(sgnum: Union[int, str]) -> List[List[float]]:
    """
    Get all unique translation parts of operations in the space group's normalizer that
    don't leave the primitive cell. Given in the primitive basis as defined in
    http://doi.org/10.1016/j.commatsci.2017.01.017

    Args:
        sgnum: space group number
    """
    _check_space_group(sgnum)
    with open(os.path.join(DATA_DIR, "possible_primitive_shifts.json")) as f:
        return json.load(f)[str(sgnum)]


def get_primitive_genpos_ops(sgnum: Union[int, str]) -> List[Dict]:
    """
    Get the matrices and columns of the space group operations in the primitive setting
    as defined in http://doi.org/10.1016/j.commatsci.2017.01.017

    Args:
        sgnum: space group number

    Returns:
        List of dictionaries, with each dictionary containing a matrix 'W' and
        translation 'w' as generally defined in the ITA, but in the primitive setting.
    """
    _check_space_group(sgnum)
    with open(os.path.join(DATA_DIR, "primitive_GENPOS_ops.json")) as f:
        return np.asarray(json.load(f)[str(sgnum)])
