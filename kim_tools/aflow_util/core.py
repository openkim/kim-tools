"""Tools for working with crystal prototypes using the AFLOW command line tool"""

import json
import logging
import os
import re
import subprocess
import sys
from curses.ascii import isalpha, isdigit
from dataclasses import dataclass
from math import acos, cos, degrees, radians, sin, sqrt
from os import PathLike
from random import random
from tempfile import NamedTemporaryFile
from typing import Any, Dict, List, Optional, Tuple, Union

import ase
import numpy as np
from ase import Atoms
from ase.cell import Cell
from numpy.typing import ArrayLike
from sympy import Symbol, linear_eq_to_matrix, matrix2numpy, parse_expr

from ..symmetry_util import (
    A_CENTERED_ORTHORHOMBIC_GROUPS,
    C_CENTERED_ORTHORHOMBIC_GROUPS,
    CENTERING_DIVISORS,
    IncorrectNumAtomsException,
    are_in_same_wyckoff_set,
    cartesian_rotation_is_in_point_group,
    get_possible_primitive_shifts,
    get_primitive_wyckoff_multiplicity,
    get_wyck_pos_xform_under_normalizer,
    space_group_numbers_are_enantiomorphic,
)

logger = logging.getLogger(__name__)
logging.basicConfig(filename="kim-tools.log", level=logging.INFO, force=True)


__author__ = ["ilia Nikiforov", "Ellad Tadmor"]
__all__ = [
    "EquivalentEqnSet",
    "EquivalentAtomSet",
    "write_tmp_poscar_from_atoms_and_run_function",
    "get_equivalent_atom_sets_from_prototype_and_atom_map",
    "IncorrectSpaceGroupException",
    "IncorrectSpeciesException",
    "InconsistentWyckoffException",
    "check_number_of_atoms",
    "split_parameter_array",
    "internal_parameter_sort_key",
    "get_stoich_reduced_list_from_prototype",
    "get_wyckoff_lists_from_prototype",
    "prototype_labels_are_equivalent",
    "get_space_group_number_from_prototype",
    "get_pearson_symbol_from_prototype",
    "get_bravais_lattice_from_prototype",
    "read_shortnames",
    "get_real_to_virtual_species_map",
    "solve_for_aflow_cell_params_from_primitive_ase_cell_params",
    "AFLOW",
]


class IncorrectSpaceGroupException(Exception):
    """
    Raised when spglib or aflow --sgdata detects a different space group than the one
    specified in the prototype label
    """


class IncorrectSpeciesException(Exception):
    """
    Raised when number or identity of species is inconsistent
    """


class InconsistentWyckoffException(Exception):
    """
    Raised when an insonsistency in Wyckoff positions is detected
    """


@dataclass
class EquivalentEqnSet:
    """
    Set of equations representing the fractional positions of equivalent atoms
    """

    species: str
    wyckoff_letter: str
    # The n free parameters associated with this Wyckoff postition, 0 <= n <= 3
    param_names: List[str]
    # m x 3 x n matrices of coefficients, where m is the multiplicity of the Wyckoff
    # position
    coeff_matrix_list: List[ArrayLike]
    # m x 3 x 1 columns of constant terms in the coordinates. This gets subtracted from
    # the RHS when solving
    const_terms_list: List[ArrayLike]


@dataclass
class EquivalentAtomSet:
    """
    Set of equivalent atoms
    """

    species: str
    wyckoff_letter: str
    frac_position_list: List[ArrayLike]  # m x 3 x 1 columns


def write_tmp_poscar_from_atoms_and_run_function(
    atoms: Atoms, function: callable, *args, **kwargs
) -> Any:
    """
    Write the Atoms file to a NamedTemporaryFile and run 'function' on it.

    Args:
        atoms:
            The atoms object that will be written to a POSCAR file and fed as the
            first argument to function
        function: A function that takes a POSCAR file as the first argument

    Returns:
        Whatever `function` returns
    """
    with NamedTemporaryFile("w+") as fp:
        atoms.write(fp, sort=True, format="vasp")
        fp.seek(0)
        return function(fp.name, *args, **kwargs)


def check_number_of_atoms(
    atoms: Atoms, prototype_label: str, primitive_cell: bool = True
) -> None:
    """
    Check if the Atoms object (which must be a conventional or primitive unit cell)
    has the correct number of atoms according to prototype_label

    Raises:
    IncorrectNumAtomsException
    """
    prototype_label_list = prototype_label.split("_")
    pearson = prototype_label_list[1]

    # get the number of atoms in conventional cell from the Pearson symbol
    num_conv_cell = 0
    for character in pearson:
        if character.isdigit():
            num_conv_cell *= 10
            num_conv_cell += int(character)

    centering = pearson[1]

    if centering == "R":
        num_conv_cell *= 3

    if not primitive_cell:
        num_lattice = 1
    else:
        num_lattice = CENTERING_DIVISORS[centering]

    # This check is probably really extraneous, but better safe than sorry
    if num_conv_cell % num_lattice != 0:
        raise IncorrectNumAtomsException(
            f"WARNING: Number of atoms in conventional cell {num_conv_cell} derived "
            f"from Pearson symbol of prototype {prototype_label} is not divisible by "
            f"the number of lattice points {num_lattice}"
        )

    num_cell = num_conv_cell / num_lattice

    if len(atoms) != num_cell:
        raise IncorrectNumAtomsException(
            f"WARNING: Number of ASE atoms {len(atoms)} does not match Pearson symbol "
            f" of prototype {prototype_label}"
        )


def split_parameter_array(
    parameter_names: List[str], list_to_split: Optional[List] = None
) -> Tuple[List, List]:
    """
    Split a list of parameters into cell and internal parameters.

    Args:
        parameter_names:
            List of AFLOW parameter names, e.g.
            `["a", "c/a", "x1", "x2", "y2", "z2"]`
            Proper AFLOW order is assumed, i.e. cell parameters
            first, then internal.
        list_to_split:
            List to split, must be same length as `parameter_names`
            If omitted, `parameter_names` itself will be split

    Returns:
        `list_to_split` (or `parameter_names` if `list_to_split` is omitted),
        split into lists corresponding to the split between cell and internal parameters

    Raises:
        AssertionError:
            If lengths are incompatible or if `parameter_names` fails an (incomplete)
            check that it is a sensible list of AFLOW parameters
    """
    if list_to_split is None:
        list_to_split = parameter_names

    assert len(list_to_split) == len(
        parameter_names
    ), "`list_to_split` must have the same length as `parameter_names`"

    in_internal_part = False

    cell_part = []
    internal_part = []

    CARTESIAN_AXES = ["x", "y", "z"]

    for name, value in zip(parameter_names, list_to_split):
        assert isinstance(
            name, str
        ), "At least one element of `parameter_names` is not a string."
        if not in_internal_part:
            if name[0] in CARTESIAN_AXES:
                in_internal_part = True
        else:
            # means we have already encountered
            # an internal coordinate in a past iteration
            assert name[0] in CARTESIAN_AXES, (
                "`parameter_names` seems to have an internal parameter "
                "followed by a non-internal one"
            )

        if in_internal_part:
            internal_part.append(value)
        else:
            cell_part.append(value)

    return cell_part, internal_part


def internal_parameter_sort_key(parameter_name: Union[Symbol, str]) -> int:
    """
    Sorting key for internal free parameters. Sort by number first, then letter
    """
    parameter_name_str = str(parameter_name)
    axis = parameter_name_str[0]
    assert (
        axis == "x" or axis == "y" or axis == "z"
    ), "Parameter name must start with x, y, or z"
    number = int(parameter_name_str[1:])
    return 1000 * number + ord(axis)


def get_equivalent_atom_sets_from_prototype_and_atom_map(
    atoms: Atoms, prototype_label: str, atom_map: List[int], sort_atoms: bool = False
) -> List[EquivalentAtomSet]:
    """
    Get a list of objects representing sets of equivalent atoms from the atom_map of an
    AFLOW comparison

    The AFLOW comparison should be between `atoms` and `atoms_rebuilt`, in that order,
    where `atoms_rebuilt` is regenerated from the prototype designation detected from
    `atoms`, and `prototype_label` is the detected prototype label

    Args:
        sort_atoms: If `atom_map` was obtained by sorting `atoms` before writing it to
        POSCAR, set this to True
    """

    if sort_atoms:
        with NamedTemporaryFile("w+") as f:
            atoms.write(f.name, format="vasp", sort=True)
            f.seek(0)
            atoms_in_correct_order = ase.io.read(f.name, format="vasp")
    else:
        atoms_in_correct_order = atoms

    # initialize return list
    equivalent_atom_set_list = []

    redetected_wyckoff_lists = get_wyckoff_lists_from_prototype(prototype_label)
    sg_num = get_space_group_number_from_prototype(prototype_label)

    index_in_atoms_rebuilt = 0
    for i_species, species_wyckoff_list in enumerate(redetected_wyckoff_lists):
        virtual_species_letter = chr(65 + i_species)
        for wyckoff_letter in species_wyckoff_list:
            equivalent_atom_set_list.append(
                EquivalentAtomSet(virtual_species_letter, wyckoff_letter, [])
            )
            for _ in range(get_primitive_wyckoff_multiplicity(sg_num, wyckoff_letter)):
                # atom_map[index_in_atoms] = index_in_atoms_rebuilt
                index_in_atoms = atom_map.index(index_in_atoms_rebuilt)
                equivalent_atom_set_list[-1].frac_position_list.append(
                    atoms_in_correct_order.get_scaled_positions()[
                        index_in_atoms
                    ].reshape(3, 1)
                )
                index_in_atoms_rebuilt += 1

    return equivalent_atom_set_list


def get_stoich_reduced_list_from_prototype(prototype_label: str) -> List[int]:
    """
    Get numerical list of stoichiometry from prototype label, i.e. "AB3_hP8..." -> [1,3]

    Args:
        prototype_label:
            AFLOW prototype label

    Returns:
        List of reduced stoichiometric numbers
    """
    stoich_reduced_formula = prototype_label.split("_")[0]
    stoich_reduced_list = []
    stoich_reduced_curr = None
    for char in stoich_reduced_formula:
        if isalpha(char):
            if stoich_reduced_curr is not None:
                if stoich_reduced_curr == 0:
                    stoich_reduced_curr = 1
                stoich_reduced_list.append(stoich_reduced_curr)
            stoich_reduced_curr = 0
        else:
            assert isdigit(char)
            # will throw an error if we haven't encountered an alphabetical letter, good
            stoich_reduced_curr *= 10
            stoich_reduced_curr += int(char)
    # write final number
    if stoich_reduced_curr == 0:
        stoich_reduced_curr = 1
    stoich_reduced_list.append(stoich_reduced_curr)
    return stoich_reduced_list


def get_wyckoff_lists_from_prototype(prototype_label: str) -> List[str]:
    """
    Expand the list of Wyckoff letters in the prototype to account for each individual
    letter instead of using numerical multipliers for repeated letters.
    e.g. A2B3C_mC48_15_aef_3f_2e -> ['aef','fff','ee']
    """
    expanded_wyckoff_letters = []
    prototype_label_split = prototype_label.split("_")
    for species_wyckoff_string in prototype_label_split[3:]:
        expanded_wyckoff_letters.append("")
        curr_wyckoff_count = 0
        for char in species_wyckoff_string:
            if isalpha(char):
                if curr_wyckoff_count == 0:
                    curr_wyckoff_count = 1
                expanded_wyckoff_letters[-1] += char * curr_wyckoff_count
                curr_wyckoff_count = 0
            else:
                assert isdigit(char)
                curr_wyckoff_count *= 10  # if it's zero, we're all good
                curr_wyckoff_count += int(char)
    return expanded_wyckoff_letters


def prototype_labels_are_equivalent(
    prototype_label_1: str,
    prototype_label_2: str,
    allow_enantiomorph: bool = False,
    allow_species_permutation: bool = False,
) -> bool:
    """
    Checks if two prototype labels are equivalent
    """
    if allow_species_permutation:
        # TODO: Possibly add this (for checking library prototype labels)
        raise NotImplementedError("Species permutations not implemented")

    if prototype_label_1 == prototype_label_2:
        return True

    # Check stoichiometry
    stoich_reduced_list_1 = get_stoich_reduced_list_from_prototype(prototype_label_1)

    if not stoich_reduced_list_1 == get_stoich_reduced_list_from_prototype(
        prototype_label_2
    ):
        logger.info(
            "Found non-matching stoichiometry in labels "
            f"{prototype_label_1} and {prototype_label_2}"
        )
        return False

    # Check Pearson symbol
    if not get_pearson_symbol_from_prototype(
        prototype_label_1
    ) == get_pearson_symbol_from_prototype(prototype_label_2):
        logger.info(
            "Found non-matching Pearson symbol in labels "
            f"{prototype_label_1} and {prototype_label_2}"
        )
        return False

    # Check space group number
    sg_num_1 = get_space_group_number_from_prototype(prototype_label_1)
    sg_num_2 = get_space_group_number_from_prototype(prototype_label_2)
    if allow_enantiomorph and not space_group_numbers_are_enantiomorphic(
        sg_num_2, sg_num_1
    ):
        logger.info(
            "Found non-matching Space group in labels "
            f"{prototype_label_1} and {prototype_label_2}"
        )
        return False
    elif sg_num_2 != sg_num_1:
        logger.info(
            "Found non-matching Space group in labels "
            f"{prototype_label_1} and {prototype_label_2}"
        )
        return False

    # OK, so far everything matches, now check the Wyckoff letters
    # Get lists of Wyckoff letters for each species,
    # e.g. A2B3C_mC48_15_aef_3f_2e -> ['aef', 'fff', 'ee']
    wyckoff_lists_1 = get_wyckoff_lists_from_prototype(prototype_label_1)
    wyckoff_lists_2 = get_wyckoff_lists_from_prototype(prototype_label_2)
    num_species = len(stoich_reduced_list_1)
    assert len(wyckoff_lists_1) == len(wyckoff_lists_2) == num_species, (
        "Somehow I got non-matching lists of Wyckoff letters, "
        "the prototype labels are probably malformed"
    )

    if sg_num_1 >= 16:
        # Theoretically, unless we are allowing species permutations, orthorhombic and
        # higher SGs should always have identical prototype labels due to minimal
        # Wyckoff enumeration. However, there are bugs in AFLOW making this untrue.

        wyck_pos_xform_list = get_wyck_pos_xform_under_normalizer(sg_num_1)
        for wyck_pos_xform in wyck_pos_xform_list:
            wyckoff_lists_match_for_each_species = True
            for i in range(num_species):
                wyckoff_list_ref = wyckoff_lists_1[i]
                wyckoff_list_test = ""
                for test_letter in wyckoff_lists_2[i]:
                    if test_letter == "A":
                        # SG47 runs out of the alphabet and uses capital A for its
                        # general position (which is unchanged under any transformation
                        # in the normalizer). However, we need "A" to be last, so we
                        # make it "{" to make it sort after all lowercase letters and
                        # replace it after
                        wyckoff_list_test += "{"
                    else:
                        test_letter_index = ord(test_letter) - 97
                        wyckoff_list_test += wyck_pos_xform[test_letter_index]
                wyckoff_list_test = "".join(sorted(wyckoff_list_test))
                wyckoff_list_test = wyckoff_list_test.replace("{", "A")
                if wyckoff_list_test != wyckoff_list_ref:
                    wyckoff_lists_match_for_each_species = False
                    break
            if wyckoff_lists_match_for_each_species:
                logger.warning(
                    f"Labels {prototype_label_1} and {prototype_label_2} were found to "
                    "be equivalent despite being non-identical. "
                    "This indicates a failure to find the lowest Wyckoff enumeration."
                )
                return True
        logger.info(
            f"Labels {prototype_label_1} and {prototype_label_2} were not found to be "
            "equivalent under any permutations allowable by the normalizer."
        )
        return False
    else:
        for wyckoff_list_1, wyckoff_list_2 in zip(wyckoff_lists_1, wyckoff_lists_2):
            for letter_1, letter_2 in zip(wyckoff_list_1, wyckoff_list_2):
                # Wyckoff sets are alphabetically contiguous in SG1-15, there is no
                # need to re-sort anything.
                # This is NOT true for all SGs (e.g. #200, Wyckoff set eh )
                if not are_in_same_wyckoff_set(letter_1, letter_2, sg_num_1):
                    logger.info(
                        f"Labels {prototype_label_1} and {prototype_label_2} have "
                        f"corresponding letters {letter_1} and {letter_2} that are not "
                        "in the same Wyckoff set"
                    )
                    return False
        logger.info(
            f"Labels {prototype_label_1} and {prototype_label_2} were found to be "
            "equivalent despite being non-identical. This is a normal occurrence for "
            "triclinic and monoclinic space groups such as this."
        )
        return True


def get_space_group_number_from_prototype(prototype_label: str) -> int:
    return int(prototype_label.split("_")[2])


def get_pearson_symbol_from_prototype(prototype_label: str) -> str:
    return prototype_label.split("_")[1]


def get_bravais_lattice_from_prototype(prototype_label: str) -> str:
    return get_pearson_symbol_from_prototype(prototype_label)[:2]


def read_shortnames() -> Dict:
    """
    This function parses "README_PROTO.TXT". It finds each line that (after stripping
    whitespace) starts with "ANRL Label". These are headers of sections of prototype
    listings. It finds the column of the word "notes". This will be the column that the
    shortnames are in. Skipping various non-prototype lines, the first column in each
    prototype line (before the ".") is the prototype, while the end of the line
    starting from the "notes" column, cleaned up to remove whitespace and
    end-of-shortname comments (i.e. "(part 3)"), is the shortname.

    Returns:
        A dictionary where the keys are the prototype strings, and the values are the
        shortnames found in the corresponding lines.
    """
    shortnames = {}
    shortname_file = "data/README_PROTO.TXT"
    notes_index = None
    with open(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), shortname_file),
        encoding="utf-8",
    ) as f:
        for line in f:
            line = line.strip()
            if line.startswith("ANRL Label"):
                try:
                    notes_index = line.index("notes")
                    continue
                except Exception:
                    print("ERROR: ANRL Label line without notes header")
                    print(line)
                    sys.exit()
            # Skip this line if it's before the first ANRL label
            if notes_index is None:
                continue
            # Skip this line if it's empty, a comment, or a divider
            if (
                line == ""
                or line == "\n"
                or line.startswith("*")
                or line.startswith("-")
                or line.startswith("ANRL")
            ):
                continue
            # Skip this line if it only has content in the first column
            # (prototype runover from previous line)
            try:
                _ = line.split(" ")[1]
            except Exception:
                continue
            # Clean up prototype (remove decorations suffix)
            prototype = line.split(" ")[0]
            if "." in prototype:
                idx = prototype.index(".")
                prototype = prototype[:idx]
            # Clean up short name
            sname = line[notes_index:]
            if "(part " in sname:
                idx = sname.index("(part")
                sname = sname[:idx]
            sname = sname.replace(", part 3", "")
            if "ICSD" in sname:
                idx = sname.index("ICSD")
                tmp = sname[idx:].split(" ")[1]
                if tmp.endswith(","):
                    sname = sname[:idx] + "ICSD " + tmp[:-1]
            if " similar to" in sname:
                idx = sname.index(" similar to")
                sname = sname[:idx]
            if " equivalent to" in sname:
                idx = sname.index(" equivalent to")
                sname = sname[:idx]
            if sname.endswith(","):
                sname = sname[:-1]
            # add prototype to shortnames dictionary
            shortnames[prototype] = sname.rstrip()
    return shortnames


def get_real_to_virtual_species_map(input: Union[List[str], Atoms]) -> Dict:
    """
    Map real species to virtual species according to (alphabetized) AFLOW convention,
    e.g. for SiC return {'C':'A', 'Si':'B'}
    """
    if isinstance(input, Atoms):
        species = sorted(list(set(input.get_chemical_symbols())))
    else:
        species = input

    real_to_virtual_species_map = {}
    for i, symbol in enumerate(species):
        real_to_virtual_species_map[symbol] = chr(65 + i)

    return real_to_virtual_species_map


def solve_for_aflow_cell_params_from_primitive_ase_cell_params(
    cellpar_prim: ArrayLike, prototype_label: str
) -> List[float]:
    """
    Get conventional cell parameters from primitive cell parameters. It is assumed that
    the primitive cell is related to the conventional cell as specified in
    10.1016/j.commatsci.2017.01.017. Equations obtained from Wolfram notebook in
    ``scripts/cell_param_solver.nb``

    Args:
        cellpar_prim:
            The 6 cell parameters of the primitive unit cell:
            [a, b, c, alpha, beta, gamma]
        prototype_label:
            The AFLOW prototype label of the crystal

    Returns:
        The cell parameters expected by AFLOW for the prototype label provided. The
        first parameter is always "a" and is given in the same units as
        ``cellpar_prim``, the others are fractional parameters in terms of "a", or
        angles in degrees. For example, if the ``prototype_label`` provided indicates a
        monoclinic crystal, this function will return the values of [a, b/a, c/a, beta]
    """
    bravais_lattice = get_bravais_lattice_from_prototype(prototype_label)

    assert len(cellpar_prim) == 6, "Got a number of cell parameters that is not 6"

    for length in cellpar_prim[0:3]:
        assert length > 0, "Got a negative cell size"
    for angle in cellpar_prim[3:]:
        assert 0 < angle < 180, "Got a cell angle outside of (0,180)"

    aprim = cellpar_prim[0]
    bprim = cellpar_prim[1]
    cprim = cellpar_prim[2]
    alphaprim = cellpar_prim[3]
    betaprim = cellpar_prim[4]
    gammaprim = cellpar_prim[5]

    if bravais_lattice == "aP":
        return [aprim, bprim / aprim, cprim / aprim, alphaprim, betaprim, gammaprim]
    elif bravais_lattice == "mP":
        return [aprim, bprim / aprim, cprim / aprim, betaprim]
    elif bravais_lattice == "oP":
        return [aprim, bprim / aprim, cprim / aprim]
    elif bravais_lattice == "tP" or bravais_lattice == "hP":
        return [aprim, cprim / aprim]
    elif bravais_lattice == "cP":
        return [aprim]
    elif bravais_lattice == "mC":
        cos_alphaprim = cos(radians(alphaprim))
        cos_gammaprim = cos(radians(gammaprim))
        a = aprim * sqrt(2 + 2 * cos_gammaprim)
        b = aprim * sqrt(2 - 2 * cos_gammaprim)
        c = cprim
        beta = degrees(acos(cos_alphaprim / sqrt((1 + cos_gammaprim) / 2)))
        return [a, b / a, c / a, beta]
    elif bravais_lattice == "oC":
        # the 'C' is colloquial, and can refer to either C or A-centering
        space_group_number = get_space_group_number_from_prototype(prototype_label)
        if space_group_number in C_CENTERED_ORTHORHOMBIC_GROUPS:
            cos_gammaprim = cos(radians(gammaprim))
            a = bprim * sqrt(2 + 2 * cos_gammaprim)
            b = bprim * sqrt(2 - 2 * cos_gammaprim)
            c = cprim
        elif space_group_number in A_CENTERED_ORTHORHOMBIC_GROUPS:
            cos_alphaprim = cos(radians(alphaprim))
            a = aprim
            b = bprim * sqrt(2 + 2 * cos_alphaprim)
            c = bprim * sqrt(2 - 2 * cos_alphaprim)
        else:
            raise IncorrectSpaceGroupException(
                f"Space group in prototype label {prototype_label} not found in lists "
                "of side-centered orthorhombic groups"
            )
        return [a, b / a, c / a]
    elif bravais_lattice == "oI":
        cos_alphaprim = cos(radians(alphaprim))
        cos_betaprim = cos(radians(betaprim))
        a = aprim * sqrt(2 + 2 * cos_alphaprim)
        b = aprim * sqrt(2 + 2 * cos_betaprim)
        # I guess the cosines must sum to a negative number!?
        # Will raise a ValueError: math domain error if not
        c = aprim * sqrt(-2 * (cos_alphaprim + cos_betaprim))
        return [a, b / a, c / a]
    elif bravais_lattice == "oF":
        aprimsq = aprim * aprim
        bprimsq = bprim * bprim
        cprimsq = cprim * cprim
        a = sqrt(2 * (-aprimsq + bprimsq + cprimsq))
        b = sqrt(2 * (aprimsq - bprimsq + cprimsq))
        c = sqrt(2 * (aprimsq + bprimsq - cprimsq))
        return [a, b / a, c / a]
    elif bravais_lattice == "tI":
        cos_alphaprim = cos(radians(alphaprim))
        a = aprim * sqrt(2 + 2 * cos_alphaprim)
        # I guess primitive alpha is always obtuse!? Will raise a
        # ValueError: math domain error if not
        c = 2 * aprim * sqrt(-cos_alphaprim)
        return [a, c / a]
    elif bravais_lattice == "hR":
        cos_alphaprim = cos(radians(alphaprim))
        a = aprim * sqrt(2 - 2 * cos_alphaprim)
        c = aprim * sqrt(3 + 6 * cos_alphaprim)
        return [a, c / a]
    elif bravais_lattice == "cF":
        return [aprim * sqrt(2)]
    elif bravais_lattice == "cI":
        return [aprim * 2 / sqrt(3)]


class AFLOW:
    """
    Class enabling access to the AFLOW executable

    Attributes:
        aflow_executable (str): Name of the AFLOW executable
        aflow_work_dir (str): Path to the work directory
        np (int):
            Number of processors to use, passed to the AFLOW executable using the
            ``--np=...`` argument

    """

    class AFLOWNotFoundException(Exception):
        """
        Raised when the AFLOW executable is not found
        """

    class ChangedSymmetryException(Exception):
        """
        Raised when an unexpected symmetry change is detected
        """

    class FailedToMatchException(Exception):
        """
        Raised when ``aflow --compare...`` fails to match
        """

    class FailedToSolveException(Exception):
        """
        Raised when solution algorithm fails
        """

    def __init__(
        self, aflow_executable: str = "aflow", aflow_work_dir: str = "", np: int = 4
    ):
        """
        Args:
            aflow_executable: Sets :attr:`aflow_executable`
            aflow_work_dir: Sets :attr:`aflow_work_dir`
            np: Sets :attr:`np`
        """
        self.aflow_executable = aflow_executable

        try:
            subprocess.check_output(["aflow", "--proto=A_cF4_225_a"])
        except Exception:
            raise self.AFLOWNotFoundException(
                "Failed to run an AFLOW test command. It is likely "
                "that the AFLOW executable was not found."
            )

        self.np = np
        if aflow_work_dir != "" and not aflow_work_dir.endswith("/"):
            self.aflow_work_dir = aflow_work_dir + "/"
        else:
            self.aflow_work_dir = aflow_work_dir

    def aflow_command(self, cmd: Union[str, List[str]], verbose=True) -> str:
        """
        Run AFLOW executable with specified arguments and return the output, possibly
        multiple times piping outputs to each other

        Args:
            cmd:
                List of arguments to pass to each AFLOW executable.
                If it's longer than 1, multiple commands will be piped to each other
            verbose: Whether to echo command to log file

        Raises:
            AFLOW.ChangedSymmetryException:
                if an ``aflow --proto=`` command complains that
                "the structure has a higher symmetry than indicated by the label"

        Returns:
            Output of the AFLOW command
        """
        if not isinstance(cmd, list):
            cmd = [cmd]

        cmd_list = [
            self.aflow_executable + " --np=" + str(self.np) + " " + cmd_inst
            for cmd_inst in cmd
        ]
        cmd_str = " | ".join(cmd_list)
        if verbose:
            logger.info(cmd_str)
        try:
            return subprocess.check_output(
                cmd_str, shell=True, stderr=subprocess.PIPE, encoding="utf-8"
            )
        except subprocess.CalledProcessError as exc:
            if "--proto=" in cmd_str and (
                "The structure has a higher symmetry than indicated by the "
                "label. The correct label and parameters for this structure are:"
            ) in str(exc.stderr):
                warn_str = (
                    "WARNING: the following command refused to write a POSCAR because "
                    f"it detected a higher symmetry: {cmd_str}. "
                    f"AFLOW error follows:\n{str(exc.stderr)}"
                )
                logger.warning(warn_str)
                raise self.ChangedSymmetryException(warn_str)
            else:
                raise exc

    def write_poscar_from_prototype(
        self,
        prototype_label: str,
        species: Optional[List[str]] = None,
        parameter_values: Optional[List[float]] = None,
        output_file: Optional[str] = None,
        verbose: bool = True,
        addtl_args: str = "",
    ) -> Optional[str]:
        """
        Run the ``aflow --proto`` command to write a POSCAR coordinate file
        corresponding to the provided AFLOW prototype designation.
        This file will have fractional coordinates.

        Args:
            prototype_label:
                An AFLOW prototype label, with or without an enumeration suffix
            species:
                List of stoichiometric species of the crystal. If this is omitted,
                the file will be written without species info
            parameter_values:
                The free parameters of the AFLOW prototype designation. If an
                enumeration suffix is not included in `prototype_label`
                and the prototype has free parameters besides `a`, this must be provided
            output_file:
                Name of the output file. If not provided,
                the output is returned as a string
            verbose: Whether to echo command to log file
            addtl_args:
                additional arguments to pass, e.g. ``--equations_only`` to get equations

        Returns:
            The output of the command or None if an `output_file` was given

        Raises:
            AFLOW.ChangedSymmetryException:
                if an ``aflow --proto=`` command complains that
                "the structure has a higher symmetry than indicated by the label"
        """
        command = " --proto=" + prototype_label
        if parameter_values:
            command += " --params=" + ",".join(
                [str(param) for param in parameter_values]
            )

        command += " " + addtl_args

        try:
            poscar_string_no_species = self.aflow_command(command, verbose=verbose)
        except self.ChangedSymmetryException as e:
            # re-raise, just indicating that this function knows about this exception
            raise e

        if species is None:
            poscar_string = poscar_string_no_species
        else:
            poscar_string = ""
            for i, line in enumerate(
                poscar_string_no_species.splitlines(keepends=True)
            ):
                poscar_string += line
                if i == 4:
                    poscar_string += " ".join(species) + "\n"
        if output_file is None:
            return poscar_string
        else:
            with open(self.aflow_work_dir + output_file, "w") as f:
                f.write(poscar_string)

    def build_atoms_from_prototype(
        self,
        prototype_label: str,
        species: List[str],
        parameter_values: Optional[List[float]] = None,
        proto_file: Optional[str] = None,
        verbose: bool = True,
    ) -> Atoms:
        """
        Build an atoms object from an AFLOW prototype designation

        Args:
            prototype_label:
                An AFLOW prototype label, with or without an enumeration suffix
            species:
                Stoichiometric species, e.g. ["Mo", "S"] corresponding to A and B
                respectively for prototype label AB2_hP6_194_c_f indicating molybdenite
            parameter_values:
                The free parameters of the AFLOW prototype designation. If an
                enumeration suffix is not included in `prototype_label`
                and the prototype has free parameters besides `a`, this must be provided
            proto_file:
                Write the POSCAR to this permanent file for debugging
                instead of a temporary file
            verbose:
                Print details in the log file

        Returns:
            Object representing unit cell of the material

        Raises:
            AFLOW.ChangedSymmetryException:
                if an ``aflow --proto=`` command complains that
                "the structure has a higher symmetry than indicated by the label"
        """
        try:
            poscar_string = self.write_poscar_from_prototype(
                prototype_label=prototype_label,
                species=species,
                parameter_values=parameter_values,
                verbose=verbose,
            )
        except self.ChangedSymmetryException as e:
            # re-raise, just indicating that this function knows about this exception
            raise e

        with (
            NamedTemporaryFile(mode="w+")
            if proto_file is None
            else open(proto_file, mode="w+")
        ) as f:
            f.write(poscar_string)
            f.seek(0)
            atoms = ase.io.read(f.name, format="vasp")
        check_number_of_atoms(atoms, prototype_label)
        atoms.wrap()
        return atoms

    def compare_materials_dir(
        self, materials_subdir: str, no_scale_volume: bool = True
    ) -> List[Dict]:
        """
        Compare a directory of materials using the aflow --compare_materials -D tool

        Args:
            materials_subdir:
                Path to the directory to compare from self.aflow_work_dir
            no_scale_volume:
                If `True`, the default behavior of allowing arbitrary scaling of
                structures before comparison is turned off

        Returns:
                Attributes of representative structures, their duplicates, and groups
                as a whole
        """
        # TODO: For efficiency, it is possible to --add_aflow_prototype_designation to
        # the representative structures. This does not help if we need duplicate
        # prototypes (for refdata), nor for library protos (as they are not ranked by
        # match like we need)
        command = " --compare_materials -D "
        command += self.aflow_work_dir + materials_subdir
        if no_scale_volume:
            command += " --no_scale_volume"
        output = self.aflow_command([command + " --screen_only --quiet --print=json"])
        res_json = json.loads(output)
        return res_json

    def get_aflow_version(self) -> str:
        """
        Run the ``aflow --version`` command to get the aflow version

        Returns:
            aflow++ executable version
        """
        command = " --version"
        output = self.aflow_command([command])
        return output.strip().split()[2]

    def compare_to_prototypes(self, input_file: str, prim: bool = True) -> List[Dict]:
        """
        Run the ``aflow --compare2prototypes`` command to compare the input structure
        to the AFLOW library of curated prototypes

        Args:
            input_file: path to the POSCAR file containing the structure to compare
            prim: whether to primitivize the structure first

        Returns:
            JSON list of dictionaries containing information about matching prototypes.
            In practice, this list should be of length zero or 1
        """
        if prim:
            command = [
                " --prim < " + self.aflow_work_dir + input_file,
                " --compare2prototypes --catalog=anrl --quiet --print=json",
            ]
        else:
            command = (
                " --compare2prototypes --catalog=anrl --quiet --print=json < "
                + self.aflow_work_dir
                + input_file
            )

        output = self.aflow_command(command)
        res_json = json.loads(output)
        return res_json

    def get_prototype_designation_from_file(
        self,
        input_file: str,
        prim: bool = True,
        force_wyckoff: bool = False,
        verbose: bool = False,
    ) -> Dict:
        """
        Run the ``aflow --prototype`` command to get the AFLOW prototype designation
            of the input structure

        Args:
            input_file: path to the POSCAR file containing the structure to analyze
            prim: whether to primitivize the structure first. Faster
            force_wyckoff:
                If the input is cif, do this to avoid re-analysis and just take the
                parameters as-is
            verbose: Whether to echo command to log file

        Returns:
            Dictionary describing the AFLOW prototype designation
            (label and parameters) of the input structure.
        """
        if prim:
            command = [
                " --prim < " + self.aflow_work_dir + input_file,
                " --prototype --print=json",
            ]
            assert not force_wyckoff, "Must specify prim=False with force_wyckoff"
        else:
            command = [
                " --prototype --print=json < " + self.aflow_work_dir + input_file
            ]

        if force_wyckoff:
            command[-1] += " --force_Wyckoff"

        output = self.aflow_command(command, verbose=verbose)
        res_json = json.loads(output)
        return res_json

    def get_prototype_designation_from_atoms(
        self, atoms: Atoms, prim: bool = True, verbose: bool = False
    ) -> Dict:
        """
        Run the ``aflow --prototype`` command to get the AFLOW prototype designation

        Args:
            atoms: atoms object to analyze
            prim: whether to primitivize the structure first
            verbose: Whether to echo command to log file

        Returns:
            Dictionary describing the AFLOW prototype designation
            (label and parameters) of the input structure.
        """
        return write_tmp_poscar_from_atoms_and_run_function(
            atoms, self.get_prototype_designation_from_file, prim=prim, verbose=verbose
        )

    def get_library_prototype_label_and_shortname_from_file(
        self, poscar_file: str, prim: bool = True, shortnames: Dict = read_shortnames()
    ) -> Tuple[Union[str, None], Union[str, None]]:
        """
        Use the aflow command line tool to determine the library prototype label for a
        structure and look up its human-readable shortname. In the case of multiple
        results, the enumeration with the smallest misfit that is in the prototypes
        list is returned. If none of the results are in the matching prototypes list,
        then the prototype with the smallest misfit is returned.

        Args:
            poscar_file:
                Path to input coordinate file
            prim: whether to primitivize the structure first
            shortnames:
                Dictionary with library prototype labels as keys and human-readable
                shortnames as values.

        Returns:
            * The library prototype label for the provided compound.
            * Shortname corresponding to this prototype
        """

        comparison_results = self.compare_to_prototypes(poscar_file, prim=prim)
        if len(comparison_results) > 1:
            # If zero results are returned it means the prototype is not in the
            # encyclopedia at all.
            # Not expecting a case where the number of results is greater than 1.
            raise RuntimeError(
                f"{comparison_results} results returned from comparison instead of "
                "zero or one as expected"
            )
        elif len(comparison_results) == 0:
            return None, None

        # Try to find the result with the smallest misfit that is in the matching
        # prototype list, otherwise return result with smallest misfit
        misfit_min_overall = 1e60
        found_overall = False
        misfit_min_inlist = 1e60
        found_inlist = False

        shortname = None
        for struct in comparison_results[0]["structures_duplicate"]:
            if struct["misfit"] < misfit_min_overall:
                misfit_min_overall = struct["misfit"]
                library_proto_overall = struct["name"]
                found_overall = True
            if struct["misfit"] < misfit_min_inlist and any(
                proto in struct["name"] for proto in shortnames
            ):
                misfit_min_inlist = struct["misfit"]
                library_proto_inlist = struct["name"]
                found_inlist = True
        if found_inlist:
            matching_library_prototype_label = library_proto_inlist
            shortname = shortnames[matching_library_prototype_label]
        elif found_overall:
            matching_library_prototype_label = library_proto_overall
        else:
            matching_library_prototype_label = None

        return matching_library_prototype_label, shortname

    def get_library_prototype_label_and_shortname_from_atoms(
        self, atoms: Atoms, prim: bool = True, shortnames: Dict = read_shortnames()
    ) -> Tuple[Union[str, None], Union[str, None]]:
        """
        Use the aflow command line tool to determine the library prototype label for a
        structure and look up its human-readable shortname. In the case of multiple
        results, the enumeration with the smallest misfit that is in the prototypes
        list is returned. If none of the results are in the matching prototypes list,
        then the prototype with the smallest misfit is returned.

        Args:
            atoms:
                Atoms object to compare
            prim: whether to primitivize the structure first
            shortnames:
                Dictionary with library prototype labels as keys and human-readable
                shortnames as values.

        Returns:
            * The library prototype label for the provided compound.
            * Shortname corresponding to this prototype
        """
        return write_tmp_poscar_from_atoms_and_run_function(
            atoms,
            self.get_library_prototype_label_and_shortname_from_file,
            prim=prim,
            shortnames=shortnames,
        )

    def _compare_poscars(self, poscar1: PathLike, poscar2: PathLike) -> Dict:
        return json.loads(
            self.aflow_command(
                [
                    f" --print=JSON --compare_materials={poscar1},{poscar2} "
                    "--screen_only --no_scale_volume --optimize_match --quiet"
                ],
                verbose=False,
            )
        )

    def _compare_Atoms(
        self,
        atoms1: Atoms,
        atoms2: Atoms,
        sort_atoms1: bool = True,
        sort_atoms2: bool = True,
    ) -> Dict:
        with NamedTemporaryFile() as f1, NamedTemporaryFile() as f2:
            atoms1.write(f1.name, "vasp", sort=sort_atoms1)
            atoms2.write(f2.name, "vasp", sort=sort_atoms2)
            f1.seek(0)
            f2.seek(0)
            compare = self._compare_poscars(f1.name, f2.name)
        return compare

    def get_basistransformation_rotation_originshift_atom_map_from_atoms(
        self,
        atoms1: Atoms,
        atoms2: Atoms,
        sort_atoms1: bool = True,
        sort_atoms2: bool = True,
    ) -> Tuple[
        Optional[ArrayLike],
        Optional[ArrayLike],
        Optional[ArrayLike],
        Optional[List[int]],
    ]:
        """
        Get operations to transform atoms2 to atoms1

        Args:
            sort_atoms1:
                Whether to sort atoms1 before comparing. If species are not
                alphabetized, this is REQUIRED. However, the `atom_map` returned will
                be w.r.t the sorted order, use with care!
            sort_atoms2: Whether to sort atoms2 before comparing.

        Returns:
            Tuple of arrays in the order:
            basis transformation, rotation, origin shift, atom_map
            atom_map[index_in_structure_1] = index_in_structure_2

        Raises:
            AFLOW.FailedToMatchException: if AFLOW fails to match the crystals
        """
        comparison_result = self._compare_Atoms(
            atoms1, atoms2, sort_atoms1, sort_atoms2
        )
        if (
            "structures_duplicate" in comparison_result[0]
            and comparison_result[0]["structures_duplicate"] != []
        ):
            return (
                np.asarray(
                    comparison_result[0]["structures_duplicate"][0][
                        "basis_transformation"
                    ]
                ),
                np.asarray(comparison_result[0]["structures_duplicate"][0]["rotation"]),
                np.asarray(
                    comparison_result[0]["structures_duplicate"][0]["origin_shift"]
                ),
                comparison_result[0]["structures_duplicate"][0]["atom_map"],
            )
        else:
            msg = "AFLOW was unable to match the provided crystals"
            logger.info(msg)
            raise self.FailedToMatchException(msg)

    def get_param_names_from_prototype(self, prototype_label: str) -> List[str]:
        """
        Get the parameter names by parsing the error message from AFLOW
        """
        try:
            self.write_poscar_from_prototype(prototype_label, parameter_values=[1.0])
            # if no exception, it's a cubic crystal with no internal free params
            return ["a"]
        except subprocess.CalledProcessError as exc:
            re_match = re.search(r"parameter_list=(.*) - \[dir", str(exc.stderr))
            if re_match is not None:
                return re_match.groups()[0].split(",")
            # if we got here, the exception didn't have the format we expected, re-raise
            raise exc

    def get_equation_sets_from_prototype(
        self, prototype_label: str
    ) -> List[EquivalentEqnSet]:
        """
        Get the symbolic equations for the fractional positions in the unit cell of an
        AFLOW prototype

        Args:
            prototype_label:
                An AFLOW prototype label, without an enumeration suffix, without
                specified atomic species

        Returns:
            List of EquivalentEqnSet objects

            Each EquivalentEqnSet contains:
                - species: The species of the atoms in this set.
                - wyckoff_letter: The Wyckoff letter corresponding to this set.
                - param_names: The names of the free parameters associated with this
                  Wyckoff position.
                - coeff_matrix_list: A list of 3 x n matrices of coefficients for the
                  free parameters.
                - const_terms_list: A list of 3 x 1 columns of constant terms in the
                  coordinates.
        """
        param_names = self.get_param_names_from_prototype(prototype_label)
        MAX_ATTEMPTS = 20
        ANGLE_NAMES = ["alpha", "beta", "gamma"]
        for i in range(MAX_ATTEMPTS):
            param_values = []
            triclinic = False
            if all([angle_name in param_names for angle_name in ANGLE_NAMES]):
                triclinic = True
                beta = 5 + random() * 165
                gamma = 5 + random() * 165
                beta_rad = radians(beta)
                gamma_rad = radians(gamma)
                max_cosalpha = abs(sin(beta_rad) * sin(gamma_rad)) + cos(
                    beta_rad
                ) * cos(gamma_rad)
                min_cosalpha = -abs(sin(beta_rad) * sin(gamma_rad)) + cos(
                    beta_rad
                ) * cos(gamma_rad)
                cosalpha = min_cosalpha + random() * (max_cosalpha - min_cosalpha)
                alpha = degrees(acos(cosalpha))
                angles_dict = {"alpha": alpha, "beta": beta, "gamma": gamma}

            for pname in param_names:
                if pname in ANGLE_NAMES:
                    if triclinic:
                        param_values.append(angles_dict[pname])
                    else:
                        param_values.append(180.0 * random())
                else:
                    param_values.append(random())

            try:
                equation_poscar = self.write_poscar_from_prototype(
                    prototype_label,
                    parameter_values=param_values,
                    addtl_args="--equations_only",
                )
                break
            except subprocess.CalledProcessError:
                if i == MAX_ATTEMPTS - 1:
                    raise RuntimeError(
                        "Random parameters failed to pass "
                        "AFLOW checks 20 times in a row"
                    )
                else:
                    pass

        # get a string with one character per Wyckoff position
        # (with possible repeated letters for positions with free params)
        wyckoff_lists = get_wyckoff_lists_from_prototype(prototype_label)
        wyckoff_joined_list = "".join(wyckoff_lists)

        coord_lines = equation_poscar.splitlines()[7:]
        coord_iter = iter(coord_lines)

        space_group_number = get_space_group_number_from_prototype(prototype_label)

        equation_sets = []

        for wyckoff_letter in wyckoff_joined_list:
            species = None  # have not seen a line yet, so don't know what species it is
            param_names = None  # same as above.
            coeff_matrix_list = []
            const_terms_list = []
            # the next n positions should be equivalent
            # corresponding to this Wyckoff position
            for _ in range(
                get_primitive_wyckoff_multiplicity(space_group_number, wyckoff_letter)
            ):
                line_split = next(coord_iter).split()
                if species is None:
                    species = line_split[3]
                elif line_split[3] != species:
                    raise InconsistentWyckoffException(
                        "Encountered different species within what I thought should be "
                        f"the lines corresponding to Wyckoff position {wyckoff_letter}"
                        f"\nEquations obtained from prototype label {prototype_label}:"
                        f"\n{equation_poscar}"
                    )
                # first, get the free parameters of this line
                curr_line_free_params = set()  # sympy.Symbol
                coordinate_expr_list = []
                for expression_string in line_split[:3]:
                    coordinate_expr = parse_expr(expression_string)
                    curr_line_free_params.update(coordinate_expr.free_symbols)
                    coordinate_expr_list.append(coordinate_expr)

                # They should all have the same number, i.e. x2, y2, z2 or x14, z14,
                # so we can just string sort them
                curr_line_free_params = list(curr_line_free_params)
                curr_line_free_params.sort(key=lambda param: str(param))

                # Each line within a Wyckoff position should
                # have the same set of free parameters
                if param_names is None:
                    param_names = curr_line_free_params
                elif param_names != curr_line_free_params:
                    raise InconsistentWyckoffException(
                        "Encountered different free params within what I thought "
                        "should be the lines corresponding to Wyckoff position "
                        f"{wyckoff_letter}\nEquations obtained from prototype label "
                        f"{prototype_label}:\n{equation_poscar}"
                    )

                # Transform to matrices and vectors
                a, b = linear_eq_to_matrix(coordinate_expr_list, param_names)

                assert a.shape == (3, len(param_names))
                assert b.shape == (3, 1)

                coeff_matrix_list.append(matrix2numpy(a, dtype=np.float64))
                const_terms_list.append(matrix2numpy(-b, dtype=np.float64))

            # Done looping over this set of equivalent positions
            equation_sets.append(
                EquivalentEqnSet(
                    species=species,
                    wyckoff_letter=wyckoff_letter,
                    param_names=[str(param_name) for param_name in param_names],
                    coeff_matrix_list=coeff_matrix_list,
                    const_terms_list=const_terms_list,
                )
            )

        # do some checks
        equation_sets_iter = iter(equation_sets)
        species = None
        for species_wyckoff_list in wyckoff_lists:
            species_must_change = True
            for _ in species_wyckoff_list:
                equation_set = next(equation_sets_iter)
                if species_must_change:
                    if equation_set.species == species:
                        raise InconsistentWyckoffException(
                            "The species in the equations obtained below are "
                            "inconsistent with the number and multiplicity of Wyckoff "
                            f"positions in prototype label {prototype_label}\n"
                            f"{equation_poscar}"
                        )
                species = equation_set.species
                species_must_change = False

        return equation_sets

    def solve_for_params_of_known_prototype(
        self,
        atoms: Atoms,
        prototype_label: str,
        max_resid: float = 1e-5,
        cell_rtol: float = 0.01,
    ) -> List[float]:
        """
        Given an Atoms object that is a primitive cell of its Bravais lattice as
        defined in doi.org/10.1016/j.commatsci.2017.01.017, and its presumed prototype
        label, solves for the free parameters of the prototype label. Raises an error
        if the solution fails (likely indicating that the Atoms object provided does
        not conform to the provided prototype label.) The Atoms object may be rotated,
        translated, and permuted, but the identity of the lattice vectors must be
        unchanged w.r.t. the crystallographic prototype. In other words, there must
        exist a permutation and translation of the fractional coordinates that enables
        them to match the equations defined by the prototype label.

        Args:
            atoms: The Atoms object to analyze
            prototype_label:
                The assumed AFLOW prototype label. `atoms` must be a primitive cell
                (with lattice vectors defined in
                doi.org/10.1016/j.commatsci.2017.01.017) conforming to the symmetry
                defined by this prototype label. Rigid body rotations, translations,
                and permutations of `atoms` relative to the AFLOW setting are allowed,
                but the identity of the lattice vectors must be unchanged.
            max_resid:
                Maximum residual allowed when attempting to match the fractional
                positions of the atoms to the crystallographic equations
            cell_rtol:
                Relative tolerance on cell lengths and angles

        Returns:
            List of free parameters that will regenerate `atoms` (up to permutations,
            rotations, and translations) when paired with `prototype_label`

        Raises:
            AFLOW.ChangedSymmetryException:
                if the symmetry of the atoms object is different from `prototype_label`
            AFLOW.FailedToMatchException:
                if AFLOW fails to match the re-generated crystal to the input crystal

        """
        # solve for cell parameters
        cell_params = solve_for_aflow_cell_params_from_primitive_ase_cell_params(
            atoms.cell.cellpar(), prototype_label
        )
        species = sorted(list(set(atoms.get_chemical_symbols())))

        # First, redetect the prototype label. We can't use this as-is because it may
        # be rotated by an operation that's within the normalizer but not
        # in the space group itself
        detected_prototype_designation = self.get_prototype_designation_from_atoms(
            atoms
        )

        prototype_label_detected = detected_prototype_designation[
            "aflow_prototype_label"
        ]

        if not prototype_labels_are_equivalent(
            prototype_label, prototype_label_detected
        ):
            msg = (
                f"Redetected prototype label {prototype_label_detected} does not match "
                f"nominal {prototype_label}."
            )
            logger.info(msg)
            raise self.ChangedSymmetryException(msg)

        # rebuild the atoms
        try:
            atoms_rebuilt = self.build_atoms_from_prototype(
                prototype_label=prototype_label_detected,
                species=species,
                parameter_values=detected_prototype_designation[
                    "aflow_prototype_params_values"
                ],
            )
        except self.ChangedSymmetryException as e:
            # re-raise, just indicating that this function knows about this exception
            raise e

        # We want the negative of the origin shift from ``atoms_rebuilt`` to ``atoms``,
        # because the origin shift is the last operation to happen, so it will be in
        # the ``atoms`` frame.

        # This function gets the transformation from its second argument to its first.
        # The origin shift is Cartesian if the POSCARs are Cartesian, which they are
        # when made from Atoms

        # Sort atoms, but do not sort atoms_rebuilt
        try:
            _, _, neg_initial_origin_shift_cart, atom_map = (
                self.get_basistransformation_rotation_originshift_atom_map_from_atoms(
                    atoms, atoms_rebuilt, sort_atoms1=True, sort_atoms2=False
                )
            )
        except self.FailedToMatchException:
            # Re-raise with a more informative error message
            msg = (
                "Execution cannot continue because AFLOW failed to match the crystal "
                "with its representation re-generated from the detected prototype "
                "designation.\nRarely, this can happen if the structure is on the edge "
                "of a symmetry increase (e.g. a BCT structure with c/a very close to 1)"
            )
            logger.error(msg)
            raise self.FailedToMatchException(msg)

        # Transpose the change of basis equation for row vectors
        initial_origin_shift_frac = (-neg_initial_origin_shift_cart) @ np.linalg.inv(
            atoms.cell
        )

        logger.info(
            "Initial shift (to SOME standard origin, not necessarily the desired one): "
            f"{-neg_initial_origin_shift_cart} (Cartesian), "
            f"{initial_origin_shift_frac} (fractional)"
        )

        position_set_list = get_equivalent_atom_sets_from_prototype_and_atom_map(
            atoms, prototype_label, atom_map, sort_atoms=True
        )

        # get equation sets
        equation_set_list = self.get_equation_sets_from_prototype(prototype_label)

        if len(position_set_list) != len(equation_set_list):
            raise InconsistentWyckoffException(
                "Number of equivalent positions detected in Atoms object "
                "did not match the number of equivalent equations given"
            )

        space_group_number = get_space_group_number_from_prototype(prototype_label)
        for prim_shift in get_possible_primitive_shifts(space_group_number):
            logger.info(
                "Additionally shifting atoms by "
                f"internal fractional translation {prim_shift}"
            )

            free_params_dict = {}
            position_set_matched_list = [False] * len(position_set_list)

            for equation_set in equation_set_list:
                # Because both equations and positions are sorted by species and
                # wyckoff letter, this should be pretty efficient
                matched_this_equation_set = False
                for i, position_set in enumerate(position_set_list):
                    if position_set_matched_list[i]:
                        continue
                    if (
                        position_set.species != equation_set.species
                    ):  # These are virtual species
                        continue
                    if not are_in_same_wyckoff_set(
                        equation_set.wyckoff_letter,
                        position_set.wyckoff_letter,
                        space_group_number,
                    ):
                        continue
                    for coeff_matrix, const_terms in zip(
                        equation_set.coeff_matrix_list, equation_set.const_terms_list
                    ):
                        for frac_position in position_set.frac_position_list:
                            # Here we use column coordinate vectors
                            frac_position_shifted = (
                                frac_position
                                + np.asarray(prim_shift).reshape(3, 1)
                                + initial_origin_shift_frac.reshape(3, 1)
                            ) % 1
                            possible_shifts = (-1, 0, 1)
                            # explore all possible shifts around zero
                            # to bring back in cell.
                            # TODO: if this is too slow (27 possibilities), write an
                            # algorithm to determine which shifts are possible
                            for shift_list in [
                                (x, y, z)
                                for x in possible_shifts
                                for y in possible_shifts
                                for z in possible_shifts
                            ]:
                                shift_array = np.asarray(shift_list).reshape(3, 1)
                                candidate_internal_param_values, resid, _, _ = (
                                    np.linalg.lstsq(
                                        coeff_matrix,
                                        frac_position_shifted
                                        - const_terms
                                        - shift_array,
                                    )
                                )
                                if len(resid) == 0 or np.max(resid) < max_resid:
                                    assert len(candidate_internal_param_values) == len(
                                        equation_set.param_names
                                    )
                                    for param_name, param_value in zip(
                                        equation_set.param_names,
                                        candidate_internal_param_values,
                                    ):
                                        assert param_name not in free_params_dict
                                        free_params_dict[param_name] = (
                                            param_value[0] % 1
                                        )  # wrap to [0,1)
                                    # should only need one to match to check off this
                                    # Wyckoff position
                                    position_set_matched_list[i] = True
                                    matched_this_equation_set = True
                                    break
                                # end loop over shifts
                            if matched_this_equation_set:
                                break
                            # end loop over positions within a position set
                        if matched_this_equation_set:
                            break
                        # end loop over equations within an equation set
                    if matched_this_equation_set:
                        break
                    # end loop over position sets
                # end loop over equation sets

            if all(position_set_matched_list):
                candidate_prototype_param_values = cell_params + [
                    free_params_dict[key]
                    for key in sorted(
                        free_params_dict.keys(), key=internal_parameter_sort_key
                    )
                ]
                # The internal shift may have taken us to an internal parameter
                # solution that represents a rotation, so we need to check
                if self.confirm_unrotated_prototype_designation(
                    atoms,
                    species,
                    prototype_label,
                    candidate_prototype_param_values,
                    cell_rtol,
                ):
                    logger.info(
                        f"Found set of parameters for prototype {prototype_label} "
                        "that is unrotated"
                    )
                    return candidate_prototype_param_values
                else:
                    logger.info(
                        f"Found set of parameters for prototype {prototype_label}, "
                        "but it was rotated relative to the original cell"
                    )
            else:
                logger.info(
                    f"Failed to solve equations for prototype {prototype_label} "
                    "on this shift attempt"
                )

        msg = (
            f"Failed to solve equations for prototype {prototype_label} "
            "on any shift attempt"
        )
        logger.info(msg)
        raise self.FailedToSolveException(msg)

    def confirm_atoms_unrotated_when_cells_aligned(
        self,
        test_atoms: Atoms,
        ref_atoms: Atoms,
        sgnum: Union[int, str],
        rtol: float = 0.01,
    ) -> bool:
        """
        Check whether `test_atoms` and `reference_atoms` are unrotated as follows:
        When the cells are in :meth:`ase.cell.Cell.standard_form`, the cells are
        identical. When both crystals are rotated to standard form (rotating the cell
        and keeping the fractional coordinates unchanged), the rotation part of the
        mapping the two crystals to each other found by AFLOW is in the point group of
        the reference crystal (using the generated crystal would give the same result).
        In other words, the crystals are identically oriented (but possibly translated)
        in reference to their lattice vectors, which in turn must be identical up to a
        rotation in reference to some Cartesian coordinate system.
        The crystals must be primitive cells as defined in
        https://doi.org/10.1016/j.commatsci.2017.01.017.

        Args:
            test_atoms:
                Primitive cell of a crystal
            ref_atoms:
                Primitive cell of a crystal
            sgnum:
                Space group number
            rtol:
                Parameter to pass to :func:`numpy.allclose` for comparing cell params.
        """
        if not np.allclose(
            ref_atoms.cell.cellpar(), test_atoms.cell.cellpar(), rtol=rtol
        ):
            logger.info(
                "Cell lengths and angles do not match.\n"
                f"Original: {ref_atoms.cell.cellpar()}\n"
                f"Regenerated: {test_atoms.cell.cellpar()}"
            )
            return False
        else:
            cell_lengths_and_angles = ref_atoms.cell.cellpar()

        test_atoms_copy = test_atoms.copy()
        ref_atoms_copy = ref_atoms.copy()

        test_atoms_copy.set_cell(
            Cell.fromcellpar(cell_lengths_and_angles), scale_atoms=True
        )
        ref_atoms_copy.set_cell(
            Cell.fromcellpar(cell_lengths_and_angles), scale_atoms=True
        )

        # the rotation below is Cartesian.
        try:
            _, cart_rot, _, _ = (
                self.get_basistransformation_rotation_originshift_atom_map_from_atoms(
                    test_atoms_copy, ref_atoms_copy
                )
            )
        except self.FailedToMatchException:
            logger.info("AFLOW failed to match the recreated crystal to reference")
            return False

        return cartesian_rotation_is_in_point_group(
            cart_rot, sgnum, test_atoms_copy.cell
        )

    def confirm_unrotated_prototype_designation(
        self,
        reference_atoms: Atoms,
        species: List[str],
        prototype_label: str,
        parameter_values: List[float],
        rtol: float = 0.01,
    ) -> bool:
        """
        Check whether the provided prototype designation recreates ``reference_atoms``
        as follows: When the cells are in :meth:`ase.cell.Cell.standard_form`, the
        cells are identical. When both crystals are rotated to standard form (rotating
        the cell and keeping the fractional coordinates unchanged), the rotation part
        of the mapping the two crystals to each other found by AFLOW is in the point
        group of the reference crystal (using the generated crystal would give the same
        result). In other words, the crystals are identically oriented (but possibly
        translated) in reference to their lattice vectors, which in turn must be
        identical up to a rotation in reference to some Cartesian coordinate system.

        Args:
            species:
                Stoichiometric species, e.g. ["Mo", "S"] corresponding to A and B
                respectively for prototype label AB2_hP6_194_c_f indicating molybdenite
            prototype_label:
                An AFLOW prototype label, without an enumeration suffix,
                without specified atomic species
            parameter_values:
                The free parameters of the AFLOW prototype designation
            rtol:
                Parameter to pass to :func:`numpy.allclose` for comparing cell params

        Returns:
            Whether or not the crystals match
        """
        test_atoms = self.build_atoms_from_prototype(
            prototype_label=prototype_label,
            species=species,
            parameter_values=parameter_values,
        )

        return self.confirm_atoms_unrotated_when_cells_aligned(
            test_atoms,
            reference_atoms,
            get_space_group_number_from_prototype(prototype_label),
            rtol,
        )
