"""
Crystal Symmetry utilities and data that are independent of AFLOW
"""

from numpy.typing import ArrayLike
from typing import Optional,List,Union
import logging
import numpy as np
from ase.atoms import Atoms


logger = logging.getLogger(__name__)
logging.basicConfig(filename='kim-tools.log',level=logging.INFO,force=True)

__all__ = [
    "KIMCrystal",
    "frac_pos_match_allow_permute_wrap",
    "frac_pos_match_allow_wrap",
    "atoms_frac_pos_match_allow_permute_wrap",
    "atoms_frac_pos_match_allow_wrap",
    "shuffle_atoms",
    "are_in_same_wyckoff_set",
    "space_group_numbers_are_enantiomorphic",
    "WYCKOFF_MULTIPLICITIES",
    "CENTERING_DIVISORS",
    "C_CENTERED_ORTHORHOMBIC_GROUPS",
    "A_CENTERED_ORTHORHOMBIC_GROUPS",
    "WYCK_POS_XFORM_UNDER_NORMALIZER",
    "SPACE_GROUPS_FOR_EACH_PEARSON",
    "get_pearson_from_space_group",
    "get_formal_pearson_from_space_group",
    "POSSIBLE_PRIMITIVE_SHIFTS"
]

C_CENTERED_ORTHORHOMBIC_GROUPS = (20,21,35,36,37,63,64,65,66,67,68)
A_CENTERED_ORTHORHOMBIC_GROUPS = (38,39,40,41)

def shuffle_atoms(
    atoms: Atoms
) -> Atoms:
    atoms_shuffled = atoms.copy()
    permute = np.random.permutation(len(atoms_shuffled))
    atoms_shuffled.set_scaled_positions([atoms.get_scaled_positions()[i] for i in permute])
    atoms_shuffled.set_chemical_symbols([atoms.get_chemical_symbols()[i] for i in permute])
    return atoms_shuffled

def _frac_pos_match_sanity_checks(
    reference_positions: ArrayLike,    
    test_positions: ArrayLike,
    reference_species: Optional[List] = None,
    test_species: Optional[List] = None,
) -> bool:
    """
    Sanity checks for comparing sets of fractional positions
    """
    if reference_species is None or test_species is None:
        if not (reference_species is None and test_species is None):
            logger.warning(
                'Refusing to compare positions when one structure has species given ' \
                'and the other does not')
            return False
        logger.info('Comparing fractional positions without species')
    else:
        if not (len(reference_positions) ==
                len(test_positions) ==
                len(reference_species) ==
                len(test_species)):
            logger.info('Atomic positions and/or species lists have different lengths between test and reference')
            return False
    
    if len(reference_positions) != len(test_positions):
        logger.info('Number of atomic positions does not match')
        return False
    return True

def frac_pos_match_allow_permute_wrap(
    reference_positions: ArrayLike,    
    test_positions: ArrayLike,
    reference_species: Optional[List] = None,
    test_species: Optional[List] = None,
) -> bool:
    """
    Check if fractional positions match allowing for permutations and PBC wrapping
    """
    if not _frac_pos_match_sanity_checks(
        reference_positions,
        test_positions,
        reference_species,
        test_species
    ):
        return False

    test_position_matched = [False]*len(test_positions)
    for i,reference_position in enumerate(reference_positions):
        for j,test_position in enumerate(test_positions):
            if test_position_matched[j]: # this position already matched. 
                continue
            position_differences = np.asarray(reference_position) - np.asarray(test_position)
            if np.allclose(position_differences,np.rint(position_differences),atol=1e-5):
                if reference_species is not None:
                    if reference_species[i] != test_species[j]:
                        logger.info(f'Reference position {i} matches test position {j} but the species do not.')
                        return False
                test_position_matched[j] = True
                break
    
    if all(test_position_matched):
        logger.info('Successfully matched the fractional positions')
        return True
    else:
        logger.info('Not all fractional positions successfully matched')
        return False

def frac_pos_match_allow_wrap(
    reference_positions: ArrayLike,    
    test_positions: ArrayLike,
    reference_species: Optional[List] = None,
    test_species: Optional[List] = None,
) -> bool:
    """
    Check if fractional positions match allowing for PBC wrapping but maintaining order
    """
    if not _frac_pos_match_sanity_checks(
        reference_positions,
        test_positions,
        reference_species,
        test_species
    ):
        return False
    
    if reference_species is not None:
        if reference_species != test_species:
            logger.info(f'Species lists do not match. Got\n{test_species}\nexpected\n{reference_species}')
            return False        
    
    for ref_pos,test_pos in zip(reference_positions, test_positions):
            position_differences = np.asarray(ref_pos) - np.asarray(test_pos)
            if not np.allclose(position_differences,np.rint(position_differences),atol=1e-5):
                logger.info(f'Failed to match positions, expected {ref_pos} got {test_pos}')
                return False
            
    logger.info('Successfully matched the fractional positions')
    return True

def atoms_frac_pos_match_allow_permute_wrap(
    reference_atoms: Atoms,
    test_atoms: Atoms
):
    return frac_pos_match_allow_permute_wrap(
        reference_atoms.get_scaled_positions(),
        test_atoms.get_scaled_positions(),
        reference_atoms.get_chemical_symbols(),
        test_atoms.get_chemical_symbols()
    )

def atoms_frac_pos_match_allow_wrap(
    reference_atoms: Atoms,
    test_atoms: Atoms
):
    return frac_pos_match_allow_wrap(
        reference_atoms.get_scaled_positions(),
        test_atoms.get_scaled_positions(),
        reference_atoms.get_chemical_symbols(),
        test_atoms.get_chemical_symbols()
    )

class KIMCrystal:
    """
    Attributes and methods of a crystal structure relevant to the KIM testing system
    
    Attributes:
        stoichiometric_species: List[str]
            List of unique species in the crystal
        prototype_label: str
            AFLOW prototype label for the crystal
        parameter_names: Optional[List[str]]
            Names of free parameters of the crystal besides 'a'. May be None if the crystal is cubic with no internal DOF.
            Should have length one less than `parameter_values_angstrom`
        parameter_values_angstrom: List[float]
            Free parameter values of the crystal. The first element in each inner list is the 'a' lattice parameter in 
            angstrom, the rest (if present) are in degrees or unitless
        library_prototype_label: Optional[str]
            AFLOW library prototype label
        short_name: Optional[List[str]]
            List of human-readable short names (e.g. "Face-Centered Cubic"), if present}
    """
    def __init__(
        prototype_label: str,
        stoichiometric_species: Optional[List[str]],
        parameter_names: Optional[List[str]],
        parameter_values_angstrom: List[float],
        library_prototype_label: Optional[str],
        short_name: Optional[List[str]],
    ):
        pass

def are_in_same_wyckoff_set(letter_1: str,letter_2: str,space_group_number: int):
    wyckoff_sets = {
    1: ["a"],
    2: ["abcdefgh", "i"],
    3: ["abcd", "e"],
    4: ["a"],
    5: ["ab", "c"],
    6: ["ab", "c"],
    7: ["a"],
    8: ["a", "b"],
    9: ["a"],
    10: ["abcdefgh", "ijkl", "mn", "o"],
    11: ["abcd", "e", "f"],
    12: ["abcd", "ef", "gh", "i", "j"],
    13: ["abcd", "ef", "g"],
    14: ["abcd", "e"],
    15: ["abcd", "e", "f"],
    16: ["abcdefgh", "ijklmnopqrst", "u"],
    17: ["abcd", "e"],
    18: ["ab", "c"],
    19: ["a"],
    20: ["ab", "c"],
    21: ["abcd", "efgh", "ij", "k", "l"],
    22: ["abcd", "efghij", "k"],
    23: ["abcd", "efghij", "k"],
    24: ["abc", "d"],
    25: ["abcd", "efgh", "i"],
    26: ["ab", "c"],
    27: ["abcd", "e"],
    28: ["ab", "c", "d"],
    29: ["a"],
    30: ["ab", "c"],
    31: ["a", "b"],
    32: ["ab", "c"],
    33: ["a"],
    34: ["ab", "c"],
    35: ["ab", "c", "de", "f"],
    36: ["a", "b"],
    37: ["ab", "c", "d"],
    38: ["ab", "c", "de", "f"],
    39: ["ab", "c", "d"],
    40: ["a", "b", "c"],
    41: ["a", "b"],
    42: ["a", "b", "cd", "e"],
    43: ["a", "b"],
    44: ["ab", "cd", "e"],
    45: ["ab", "c"],
    46: ["a", "b", "c"],
    47: ["abcdefgh", "ijklmnopqrst", "uvwxyz", "A"],
    48: ["abcd", "ef", "ghijkl", "m"],
    49: ["abcd", "efgh", "ijkl", "mnop", "q", "r"],
    50: ["abcd", "ef", "ghij", "kl", "m"],
    51: ["abcd", "ef", "gh", "ij", "k", "l"],
    52: ["ab", "c", "d", "e"],
    53: ["abcd", "ef", "g", "h", "i"],
    54: ["ab", "c", "de", "f"],
    55: ["abcd", "ef", "gh", "i"],
    56: ["ab", "cd", "e"],
    57: ["ab", "c", "d", "e"],
    58: ["abcd", "ef", "g", "h"],
    59: ["ab", "cd", "ef", "g"],
    60: ["ab", "c", "d"],
    61: ["ab", "c"],
    62: ["ab", "c", "d"],
    63: ["ab", "c", "d", "e", "f", "g", "h"],
    64: ["ab", "c", "d", "e", "f", "g"],
    65: ["abcd", "ef", "ghij", "kl", "m", "no", "pq", "r"],
    66: ["ab", "cd", "ef", "gh", "ij", "k", "l", "m"],
    67: ["ab", "cdef", "g", "hijk", "l", "mn", "o"],
    68: ["ab", "cd", "ef", "g", "h", "i"],
    69: ["ab", "cde", "f", "ghi", "jkl", "mno", "p"],
    70: ["ab", "cd", "efg", "h"],
    71: ["abcd", "efghij", "k", "lmn", "o"],
    72: ["ab", "cd", "e", "fg", "hi", "j", "k"],
    73: ["ab", "cde", "f"],
    74: ["abcd", "e", "fg", "hi", "j"],
    75: ["ab", "c", "d"],
    76: ["a"],
    77: ["ab", "c", "d"],
    78: ["a"],
    79: ["a", "b", "c"],
    80: ["a", "b"],
    81: ["abcd", "ef", "g", "h"],
    82: ["abcd", "ef", "g"],
    83: ["abcd", "ef", "gh", "i", "jk", "l"],
    84: ["ab", "cd", "ef", "gh", "i", "j", "k"],
    85: ["ab", "c", "de", "f", "g"],
    86: ["ab", "cd", "e", "f", "g"],
    87: ["ab", "c", "d", "e", "f", "g", "h", "i"],
    88: ["ab", "cd", "e", "f"],
    89: ["abcd", "ef", "gh", "i", "jk", "lmno", "p"],
    90: ["ab", "c", "d", "ef", "g"],
    91: ["ab", "c", "d"],
    92: ["a", "b"],
    93: ["ab", "cd", "ef", "gh", "i", "jklm", "no", "p"],
    94: ["ab", "c", "d", "ef", "g"],
    95: ["ab", "c", "d"],
    96: ["a", "b"],
    97: ["ab", "c", "d", "e", "f", "g", "hi", "j", "k"],
    98: ["ab", "c", "de", "f", "g"],
    99: ["ab", "c", "d", "ef", "g"],
    100: ["a", "b", "c", "d"],
    101: ["ab", "c", "d", "e"],
    102: ["a", "b", "c", "d"],
    103: ["ab", "c", "d"],
    104: ["a", "b", "c"],
    105: ["ab", "c", "de", "f"],
    106: ["a", "b", "c"],
    107: ["a", "b", "c", "d", "e"],
    108: ["a", "b", "c", "d"],
    109: ["a", "b", "c"],
    110: ["a", "b"],
    111: ["abcd", "ef", "gh", "ijkl", "m", "n", "o"],
    112: ["ac", "ef", "bd", "ghij", "kl", "m", "n"],
    113: ["ab", "c", "d", "e", "f"],
    114: ["ab", "c", "d", "e"],
    115: ["abcd", "ef", "g", "hi", "jk", "l"],
    116: ["ab", "cd", "ef", "gh", "i", "j"],
    117: ["ab", "cd", "e", "f", "gh", "i"],
    118: ["ab", "cd", "e", "fg", "h", "i"],
    119: ["abcd", "ef", "gh", "i", "j"],
    120: ["ad", "bc", "eh", "fg", "i"],
    121: ["ab", "c", "d", "e", "fg", "h", "i", "j"],
    122: ["ab", "c", "d", "e"],
    123: ["abcd", "ef", "gh", "i", "jk", "lmno", "pq", "r", "st", "u"],
    124: ["ac", "bd", "e", "f", "gh", "i", "j", "kl", "m", "n"],
    125: ["ab", "cd", "ef", "g", "h", "ij", "kl", "m", "n"],
    126: ["ab", "c", "d", "e", "f", "g", "h", "ij", "k"],
    127: ["ab", "cd", "e", "f", "gh", "ij", "k", "l"],
    128: ["ab", "c", "d", "e", "f", "g", "h", "i"],
    129: ["ab", "c", "de", "f", "gh", "i", "j", "k"],
    130: ["a", "b", "c", "d", "e", "f", "g"],
    131: ["ab", "cd", "ef", "gh", "i", "jklm", "n", "op", "q", "r"],
    132: ["ac", "bd", "e", "f", "gh", "k", "ij", "lm", "n", "o", "p"],
    133: ["a", "b", "c", "d", "e", "f", "g", "hi", "j", "k"],
    134: ["ab", "c", "d", "ef", "g", "h", "ij", "kl", "m", "n"],
    135: ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
    136: ["ab", "c", "d", "e", "fg", "h", "i", "j", "k"],
    137: ["ab", "c", "d", "e", "f", "g", "h"],
    138: ["a", "b", "cd", "e", "f", "gh", "i", "j"],
    139: ["ab", "c", "d", "e", "f", "g", "h", "ij", "k", "l", "m", "n", "o"],
    140: ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"],
    141: ["ab", "cd", "e", "f", "g", "h", "i"],
    142: ["a", "b", "c", "d", "e", "f", "g"],
    143: ["abc", "d"],
    144: ["a"],
    145: ["a"],
    146: ["a", "b"],
    147: ["ab", "c", "d", "ef", "g"],
    148: ["ab", "c", "de", "f"],
    149: ["abcdef", "ghi", "jk", "l"],
    150: ["ab", "c", "d", "ef", "g"],
    151: ["ab", "c"],
    152: ["ab", "c"],
    153: ["ab", "c"],
    154: ["ab", "c"],
    155: ["ab", "c", "de", "f"],
    156: ["abc", "d", "e"],
    157: ["a", "b", "c", "d"],
    158: ["abc", "d"],
    159: ["a", "b", "c"],
    160: ["a", "b", "c"],
    161: ["a", "b"],
    162: ["ab", "cd", "e", "fg", "h", "ij", "k", "l"],
    163: ["a", "b", "cd", "e", "f", "g", "h", "i"],
    164: ["ab", "c", "d", "ef", "gh", "i", "j"],
    165: ["a", "b", "c", "d", "e", "f", "g"],
    166: ["ab", "c", "de", "fg", "h", "i"],
    167: ["a", "b", "c", "d", "e", "f"],
    168: ["a", "b", "c", "d"],
    169: ["a"],
    170: ["a"],
    171: ["a", "b", "c"],
    172: ["a", "b", "c"],
    173: ["a", "b", "c"],
    174: ["abcdef", "ghi", "jk", "l"],
    175: ["ab", "cd", "e", "fg", "h", "i", "jk", "l"],
    176: ["a", "b", "cd", "e", "f", "g", "h", "i"],
    177: ["ab", "cd", "e", "fg", "h", "i", "jk", "lm", "n"],
    178: ["a", "b", "c"],
    179: ["a", "b", "c"],
    180: ["ab", "cd", "e", "f", "gh", "ij", "k"],
    181: ["ab", "cd", "e", "f", "gh", "ij", "k"],
    182: ["a", "b", "cd", "e", "f", "g", "h", "i"],
    183: ["a", "b", "c", "d", "e", "f"],
    184: ["a", "b", "c", "d"],
    185: ["a", "b", "c", "d"],
    186: ["a", "b", "c", "d"],
    187: ["abcdef", "ghi", "jk", "lm", "n", "o"],
    188: ["ace", "bdf", "ghi", "j", "k", "l"],
    189: ["ab", "cd", "e", "fg", "h", "i", "jk", "l"],
    190: ["a", "b", "cd", "e", "f", "g", "h", "i"],
    191: ["ab", "cd", "e", "fg", "h", "i", "jk", "lm", "n", "o", "pq", "r"],
    192: ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"],
    193: ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
    194: ["a", "b", "cd", "e", "f", "g", "h", "i", "j", "k", "l"],
    195: ["ab", "cd", "e", "fi", "gh", "j"],
    196: ["abcd", "e", "fg", "h"],
    197: ["a", "b", "c", "d", "e", "f"],
    198: ["a", "b"],
    199: ["a", "b", "c"],
    200: ["ab", "cd", "eh", "fg", "i", "jk", "l"],
    201: ["a", "bc", "d", "e", "f", "g", "h"],
    202: ["ab", "c", "d", "e", "f", "g", "h", "i"],
    203: ["ab", "cd", "e", "f", "g"],
    204: ["a", "b", "c", "d", "e", "f", "g", "h"],
    205: ["ab", "c", "d"],
    206: ["ab", "c", "d", "e"],
    207: ["ab", "cd", "ef", "g", "h", "ij", "k"],
    208: ["a", "bc", "d", "ef", "g", "h", "ij", "kl", "m"],
    209: ["ab", "c", "d", "e", "f", "gh", "i", "j"],
    210: ["ab", "cd", "e", "f", "g", "h"],
    211: ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
    212: ["ab", "c", "d", "e"],
    213: ["ab", "c", "d", "e"],
    214: ["ab", "cd", "e", "f", "gh", "i"],
    215: ["ab", "cd", "e", "fg", "h", "i", "j"],
    216: ["abcd", "e", "fg", "h", "i"],
    217: ["a", "b", "c", "d", "e", "f", "g", "h"],
    218: ["a", "b", "cd", "e", "f", "gh", "i"],
    219: ["ab", "cd", "e", "fg", "h"],
    220: ["ab", "c", "d", "e"],
    221: ["ab", "cd", "ef", "g", "h", "ij", "kl", "m", "n"],
    222: ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
    223: ["a", "b", "cd", "e", "f", "gh", "i", "j", "k", "l"],
    224: ["a", "bc", "d", "e", "f", "g", "h", "ij", "k", "l"],
    225: ["ab", "c", "d", "e", "f", "g", "hi", "j", "k", "l"],
    226: ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
    227: ["ab", "cd", "e", "f", "g", "h", "i"],
    228: ["a", "b", "c", "d", "e", "f", "g", "h"],
    229: ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
    230: ["a", "b", "c", "d", "e", "f", "g", "h"],
    }
    for wyckoff_set in wyckoff_sets[space_group_number]:
        if letter_1 in wyckoff_set:
            if letter_2 in wyckoff_set:
                return True
            else:
                return False

def space_group_numbers_are_enantiomorphic(sg_1: int, sg_2: int) -> bool:
    if sg_1 == sg_2:
        return True
    else:
        enantiomorph_conversion = {78:76,95:91,96:92,145:144,153:151,154:152,170:169,172:171,179:178,181:180,213:212}
        enantiomorph_conversion_2 = {v: k for k, v in enantiomorph_conversion.items()}
        enantiomorph_conversion.update(enantiomorph_conversion_2)
        if enantiomorph_conversion[sg_1] == sg_2:
            return True
        else:
            return False

WYCKOFF_MULTIPLICITIES = {
    1: {
        "a": 1,
    },
    2: {
        "i": 2,
        "h": 1,
        "g": 1,
        "f": 1,
        "e": 1,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    3: {
        "e": 2,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    4: {
        "a": 2,
    },
    5: {
        "c": 4,
        "b": 2,
        "a": 2,
    },
    6: {
        "c": 2,
        "b": 1,
        "a": 1,
    },
    7: {
        "a": 2,
    },
    8: {
        "b": 4,
        "a": 2,
    },
    9: {
        "a": 4,
    },
    10: {
        "o": 4,
        "n": 2,
        "m": 2,
        "l": 2,
        "k": 2,
        "j": 2,
        "i": 2,
        "h": 1,
        "g": 1,
        "f": 1,
        "e": 1,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    11: {
        "f": 4,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    12: {
        "j": 8,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    13: {
        "g": 4,
        "f": 2,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    14: {
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    15: {
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    16: {
        "u": 4,
        "t": 2,
        "s": 2,
        "r": 2,
        "q": 2,
        "p": 2,
        "o": 2,
        "n": 2,
        "m": 2,
        "l": 2,
        "k": 2,
        "j": 2,
        "i": 2,
        "h": 1,
        "g": 1,
        "f": 1,
        "e": 1,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    17: {
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    18: {
        "c": 4,
        "b": 2,
        "a": 2,
    },
    19: {
        "a": 4,
    },
    20: {
        "c": 8,
        "b": 4,
        "a": 4,
    },
    21: {
        "l": 8,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    22: {
        "k": 16,
        "j": 8,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    23: {
        "k": 8,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    24: {
        "d": 8,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    25: {
        "i": 4,
        "h": 2,
        "g": 2,
        "f": 2,
        "e": 2,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    26: {
        "c": 4,
        "b": 2,
        "a": 2,
    },
    27: {
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    28: {
        "d": 4,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    29: {
        "a": 4,
    },
    30: {
        "c": 4,
        "b": 2,
        "a": 2,
    },
    31: {
        "b": 4,
        "a": 2,
    },
    32: {
        "c": 4,
        "b": 2,
        "a": 2,
    },
    33: {
        "a": 4,
    },
    34: {
        "c": 4,
        "b": 2,
        "a": 2,
    },
    35: {
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    36: {
        "b": 8,
        "a": 4,
    },
    37: {
        "d": 8,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    38: {
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    39: {
        "d": 8,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    40: {
        "c": 8,
        "b": 4,
        "a": 4,
    },
    41: {
        "b": 8,
        "a": 4,
    },
    42: {
        "e": 16,
        "d": 8,
        "c": 8,
        "b": 8,
        "a": 4,
    },
    43: {
        "b": 16,
        "a": 8,
    },
    44: {
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    45: {
        "c": 8,
        "b": 4,
        "a": 4,
    },
    46: {
        "c": 8,
        "b": 4,
        "a": 4,
    },
    47: {
        "A": 8,
        "z": 4,
        "y": 4,
        "x": 4,
        "w": 4,
        "v": 4,
        "u": 4,
        "t": 2,
        "s": 2,
        "r": 2,
        "q": 2,
        "p": 2,
        "o": 2,
        "n": 2,
        "m": 2,
        "l": 2,
        "k": 2,
        "j": 2,
        "i": 2,
        "h": 1,
        "g": 1,
        "f": 1,
        "e": 1,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    48: {
        "m": 8,
        "l": 4,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    49: {
        "r": 8,
        "q": 4,
        "p": 4,
        "o": 4,
        "n": 4,
        "m": 4,
        "l": 4,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 2,
        "g": 2,
        "f": 2,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    50: {
        "m": 8,
        "l": 4,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    51: {
        "l": 8,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 2,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    52: {
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    53: {
        "i": 8,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    54: {
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    55: {
        "i": 8,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    56: {
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    57: {
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    58: {
        "h": 8,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    59: {
        "g": 8,
        "f": 4,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    60: {
        "d": 8,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    61: {
        "c": 8,
        "b": 4,
        "a": 4,
    },
    62: {
        "d": 8,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    63: {
        "h": 16,
        "g": 8,
        "f": 8,
        "e": 8,
        "d": 8,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    64: {
        "g": 16,
        "f": 8,
        "e": 8,
        "d": 8,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    65: {
        "r": 16,
        "q": 8,
        "p": 8,
        "o": 8,
        "n": 8,
        "m": 8,
        "l": 4,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    66: {
        "m": 16,
        "l": 8,
        "k": 8,
        "j": 8,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 4,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    67: {
        "o": 16,
        "n": 8,
        "m": 8,
        "l": 8,
        "k": 8,
        "j": 8,
        "i": 8,
        "h": 8,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    68: {
        "i": 16,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 8,
        "d": 8,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    69: {
        "p": 32,
        "o": 16,
        "n": 16,
        "m": 16,
        "l": 16,
        "k": 16,
        "j": 16,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 8,
        "d": 8,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    70: {
        "h": 32,
        "g": 16,
        "f": 16,
        "e": 16,
        "d": 16,
        "c": 16,
        "b": 8,
        "a": 8,
    },
    71: {
        "o": 16,
        "n": 8,
        "m": 8,
        "l": 8,
        "k": 8,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    72: {
        "k": 16,
        "j": 8,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    73: {
        "f": 16,
        "e": 8,
        "d": 8,
        "c": 8,
        "b": 8,
        "a": 8,
    },
    74: {
        "j": 16,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    75: {
        "d": 4,
        "c": 2,
        "b": 1,
        "a": 1,
    },
    76: {
        "a": 4,
    },
    77: {
        "d": 4,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    78: {
        "a": 4,
    },
    79: {
        "c": 8,
        "b": 4,
        "a": 2,
    },
    80: {
        "b": 8,
        "a": 4,
    },
    81: {
        "h": 4,
        "g": 2,
        "f": 2,
        "e": 2,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    82: {
        "g": 8,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    83: {
        "l": 8,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 2,
        "g": 2,
        "f": 2,
        "e": 2,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    84: {
        "k": 8,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 2,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    85: {
        "g": 8,
        "f": 4,
        "e": 4,
        "d": 4,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    86: {
        "g": 8,
        "f": 4,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    87: {
        "i": 16,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    88: {
        "f": 16,
        "e": 8,
        "d": 8,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    89: {
        "p": 8,
        "o": 4,
        "n": 4,
        "m": 4,
        "l": 4,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 2,
        "g": 2,
        "f": 2,
        "e": 2,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    90: {
        "g": 8,
        "f": 4,
        "e": 4,
        "d": 4,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    91: {
        "d": 8,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    92: {
        "b": 8,
        "a": 4,
    },
    93: {
        "p": 8,
        "o": 4,
        "n": 4,
        "m": 4,
        "l": 4,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 2,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    94: {
        "g": 8,
        "f": 4,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    95: {
        "d": 8,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    96: {
        "b": 8,
        "a": 4,
    },
    97: {
        "k": 16,
        "j": 8,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    98: {
        "g": 16,
        "f": 8,
        "e": 8,
        "d": 8,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    99: {
        "g": 8,
        "f": 4,
        "e": 4,
        "d": 4,
        "c": 2,
        "b": 1,
        "a": 1,
    },
    100: {
        "d": 8,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    101: {
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    102: {
        "d": 8,
        "c": 4,
        "b": 4,
        "a": 2,
    },
    103: {
        "d": 8,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    104: {
        "c": 8,
        "b": 4,
        "a": 2,
    },
    105: {
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    106: {
        "c": 8,
        "b": 4,
        "a": 4,
    },
    107: {
        "e": 16,
        "d": 8,
        "c": 8,
        "b": 4,
        "a": 2,
    },
    108: {
        "d": 16,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    109: {
        "c": 16,
        "b": 8,
        "a": 4,
    },
    110: {
        "b": 16,
        "a": 8,
    },
    111: {
        "o": 8,
        "n": 4,
        "m": 4,
        "l": 4,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 2,
        "g": 2,
        "f": 2,
        "e": 2,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    112: {
        "n": 8,
        "m": 4,
        "l": 4,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 2,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    113: {
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    114: {
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    115: {
        "l": 8,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 2,
        "f": 2,
        "e": 2,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    116: {
        "j": 8,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    117: {
        "i": 8,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    118: {
        "i": 8,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    119: {
        "j": 16,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    120: {
        "i": 16,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    121: {
        "j": 16,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    122: {
        "e": 16,
        "d": 8,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    123: {
        "u": 16,
        "t": 8,
        "s": 8,
        "r": 8,
        "q": 8,
        "p": 8,
        "o": 4,
        "n": 4,
        "m": 4,
        "l": 4,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 2,
        "g": 2,
        "f": 2,
        "e": 2,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    124: {
        "n": 16,
        "m": 8,
        "l": 8,
        "k": 8,
        "j": 8,
        "i": 8,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    125: {
        "n": 16,
        "m": 8,
        "l": 8,
        "k": 8,
        "j": 8,
        "i": 8,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    126: {
        "k": 16,
        "j": 8,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    127: {
        "l": 16,
        "k": 8,
        "j": 8,
        "i": 8,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    128: {
        "i": 16,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    129: {
        "k": 16,
        "j": 8,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 4,
        "e": 4,
        "d": 4,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    130: {
        "g": 16,
        "f": 8,
        "e": 8,
        "d": 8,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    131: {
        "r": 16,
        "q": 8,
        "p": 8,
        "o": 8,
        "n": 8,
        "m": 4,
        "l": 4,
        "k": 4,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 2,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    132: {
        "p": 16,
        "o": 8,
        "n": 8,
        "m": 8,
        "l": 8,
        "k": 8,
        "j": 4,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    133: {
        "k": 16,
        "j": 8,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    134: {
        "n": 16,
        "m": 8,
        "l": 8,
        "k": 8,
        "j": 8,
        "i": 8,
        "h": 8,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    135: {
        "i": 16,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    136: {
        "k": 16,
        "j": 8,
        "i": 8,
        "h": 8,
        "g": 4,
        "f": 4,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    137: {
        "h": 16,
        "g": 8,
        "f": 8,
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    138: {
        "j": 16,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    139: {
        "o": 32,
        "n": 16,
        "m": 16,
        "l": 16,
        "k": 16,
        "j": 8,
        "i": 8,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    140: {
        "m": 32,
        "l": 16,
        "k": 16,
        "j": 16,
        "i": 16,
        "h": 8,
        "g": 8,
        "f": 8,
        "e": 8,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    141: {
        "i": 32,
        "h": 16,
        "g": 16,
        "f": 16,
        "e": 8,
        "d": 8,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    142: {
        "g": 32,
        "f": 16,
        "e": 16,
        "d": 16,
        "c": 16,
        "b": 8,
        "a": 8,
    },
    143: {
        "d": 3,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    144: {
        "a": 3,
    },
    145: {
        "a": 3,
    },
    146: {
        "b": 9,
        "a": 3,
    },
    147: {
        "g": 6,
        "f": 3,
        "e": 3,
        "d": 2,
        "c": 2,
        "b": 1,
        "a": 1,
    },
    148: {
        "f": 18,
        "e": 9,
        "d": 9,
        "c": 6,
        "b": 3,
        "a": 3,
    },
    149: {
        "l": 6,
        "k": 3,
        "j": 3,
        "i": 2,
        "h": 2,
        "g": 2,
        "f": 1,
        "e": 1,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    150: {
        "g": 6,
        "f": 3,
        "e": 3,
        "d": 2,
        "c": 2,
        "b": 1,
        "a": 1,
    },
    151: {
        "c": 6,
        "b": 3,
        "a": 3,
    },
    152: {
        "c": 6,
        "b": 3,
        "a": 3,
    },
    153: {
        "c": 6,
        "b": 3,
        "a": 3,
    },
    154: {
        "c": 6,
        "b": 3,
        "a": 3,
    },
    155: {
        "f": 18,
        "e": 9,
        "d": 9,
        "c": 6,
        "b": 3,
        "a": 3,
    },
    156: {
        "e": 6,
        "d": 3,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    157: {
        "d": 6,
        "c": 3,
        "b": 2,
        "a": 1,
    },
    158: {
        "d": 6,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    159: {
        "c": 6,
        "b": 2,
        "a": 2,
    },
    160: {
        "c": 18,
        "b": 9,
        "a": 3,
    },
    161: {
        "b": 18,
        "a": 6,
    },
    162: {
        "l": 12,
        "k": 6,
        "j": 6,
        "i": 6,
        "h": 4,
        "g": 3,
        "f": 3,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 1,
        "a": 1,
    },
    163: {
        "i": 12,
        "h": 6,
        "g": 6,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    164: {
        "j": 12,
        "i": 6,
        "h": 6,
        "g": 6,
        "f": 3,
        "e": 3,
        "d": 2,
        "c": 2,
        "b": 1,
        "a": 1,
    },
    165: {
        "g": 12,
        "f": 6,
        "e": 6,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    166: {
        "i": 36,
        "h": 18,
        "g": 18,
        "f": 18,
        "e": 9,
        "d": 9,
        "c": 6,
        "b": 3,
        "a": 3,
    },
    167: {
        "f": 36,
        "e": 18,
        "d": 18,
        "c": 12,
        "b": 6,
        "a": 6,
    },
    168: {
        "d": 6,
        "c": 3,
        "b": 2,
        "a": 1,
    },
    169: {
        "a": 6,
    },
    170: {
        "a": 6,
    },
    171: {
        "c": 6,
        "b": 3,
        "a": 3,
    },
    172: {
        "c": 6,
        "b": 3,
        "a": 3,
    },
    173: {
        "c": 6,
        "b": 2,
        "a": 2,
    },
    174: {
        "l": 6,
        "k": 3,
        "j": 3,
        "i": 2,
        "h": 2,
        "g": 2,
        "f": 1,
        "e": 1,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    175: {
        "l": 12,
        "k": 6,
        "j": 6,
        "i": 6,
        "h": 4,
        "g": 3,
        "f": 3,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 1,
        "a": 1,
    },
    176: {
        "i": 12,
        "h": 6,
        "g": 6,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    177: {
        "n": 12,
        "m": 6,
        "l": 6,
        "k": 6,
        "j": 6,
        "i": 6,
        "h": 4,
        "g": 3,
        "f": 3,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 1,
        "a": 1,
    },
    178: {
        "c": 12,
        "b": 6,
        "a": 6,
    },
    179: {
        "c": 12,
        "b": 6,
        "a": 6,
    },
    180: {
        "k": 12,
        "j": 6,
        "i": 6,
        "h": 6,
        "g": 6,
        "f": 6,
        "e": 6,
        "d": 3,
        "c": 3,
        "b": 3,
        "a": 3,
    },
    181: {
        "k": 12,
        "j": 6,
        "i": 6,
        "h": 6,
        "g": 6,
        "f": 6,
        "e": 6,
        "d": 3,
        "c": 3,
        "b": 3,
        "a": 3,
    },
    182: {
        "i": 12,
        "h": 6,
        "g": 6,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    183: {
        "f": 12,
        "e": 6,
        "d": 6,
        "c": 3,
        "b": 2,
        "a": 1,
    },
    184: {
        "d": 12,
        "c": 6,
        "b": 4,
        "a": 2,
    },
    185: {
        "d": 12,
        "c": 6,
        "b": 4,
        "a": 2,
    },
    186: {
        "d": 12,
        "c": 6,
        "b": 2,
        "a": 2,
    },
    187: {
        "o": 12,
        "n": 6,
        "m": 6,
        "l": 6,
        "k": 3,
        "j": 3,
        "i": 2,
        "h": 2,
        "g": 2,
        "f": 1,
        "e": 1,
        "d": 1,
        "c": 1,
        "b": 1,
        "a": 1,
    },
    188: {
        "l": 12,
        "k": 6,
        "j": 6,
        "i": 4,
        "h": 4,
        "g": 4,
        "f": 2,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    189: {
        "l": 12,
        "k": 6,
        "j": 6,
        "i": 6,
        "h": 4,
        "g": 3,
        "f": 3,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 1,
        "a": 1,
    },
    190: {
        "i": 12,
        "h": 6,
        "g": 6,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    191: {
        "r": 24,
        "q": 12,
        "p": 12,
        "o": 12,
        "n": 12,
        "m": 6,
        "l": 6,
        "k": 6,
        "j": 6,
        "i": 6,
        "h": 4,
        "g": 3,
        "f": 3,
        "e": 2,
        "d": 2,
        "c": 2,
        "b": 1,
        "a": 1,
    },
    192: {
        "m": 24,
        "l": 12,
        "k": 12,
        "j": 12,
        "i": 12,
        "h": 8,
        "g": 6,
        "f": 6,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    193: {
        "l": 24,
        "k": 12,
        "j": 12,
        "i": 12,
        "h": 8,
        "g": 6,
        "f": 6,
        "e": 4,
        "d": 4,
        "c": 4,
        "b": 2,
        "a": 2,
    },
    194: {
        "l": 24,
        "k": 12,
        "j": 12,
        "i": 12,
        "h": 6,
        "g": 6,
        "f": 4,
        "e": 4,
        "d": 2,
        "c": 2,
        "b": 2,
        "a": 2,
    },
    195: {
        "j": 12,
        "i": 6,
        "h": 6,
        "g": 6,
        "f": 6,
        "e": 4,
        "d": 3,
        "c": 3,
        "b": 1,
        "a": 1,
    },
    196: {
        "h": 48,
        "g": 24,
        "f": 24,
        "e": 16,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    197: {
        "f": 24,
        "e": 12,
        "d": 12,
        "c": 8,
        "b": 6,
        "a": 2,
    },
    198: {
        "b": 12,
        "a": 4,
    },
    199: {
        "c": 24,
        "b": 12,
        "a": 8,
    },
    200: {
        "l": 24,
        "k": 12,
        "j": 12,
        "i": 8,
        "h": 6,
        "g": 6,
        "f": 6,
        "e": 6,
        "d": 3,
        "c": 3,
        "b": 1,
        "a": 1,
    },
    201: {
        "h": 24,
        "g": 12,
        "f": 12,
        "e": 8,
        "d": 6,
        "c": 4,
        "b": 4,
        "a": 2,
    },
    202: {
        "i": 96,
        "h": 48,
        "g": 48,
        "f": 32,
        "e": 24,
        "d": 24,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    203: {
        "g": 96,
        "f": 48,
        "e": 32,
        "d": 16,
        "c": 16,
        "b": 8,
        "a": 8,
    },
    204: {
        "h": 48,
        "g": 24,
        "f": 16,
        "e": 12,
        "d": 12,
        "c": 8,
        "b": 6,
        "a": 2,
    },
    205: {
        "d": 24,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    206: {
        "e": 48,
        "d": 24,
        "c": 16,
        "b": 8,
        "a": 8,
    },
    207: {
        "k": 24,
        "j": 12,
        "i": 12,
        "h": 12,
        "g": 8,
        "f": 6,
        "e": 6,
        "d": 3,
        "c": 3,
        "b": 1,
        "a": 1,
    },
    208: {
        "m": 24,
        "l": 12,
        "k": 12,
        "j": 12,
        "i": 12,
        "h": 12,
        "g": 8,
        "f": 6,
        "e": 6,
        "d": 6,
        "c": 4,
        "b": 4,
        "a": 2,
    },
    209: {
        "j": 96,
        "i": 48,
        "h": 48,
        "g": 48,
        "f": 32,
        "e": 24,
        "d": 24,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    210: {
        "h": 96,
        "g": 48,
        "f": 48,
        "e": 32,
        "d": 16,
        "c": 16,
        "b": 8,
        "a": 8,
    },
    211: {
        "j": 48,
        "i": 24,
        "h": 24,
        "g": 24,
        "f": 16,
        "e": 12,
        "d": 12,
        "c": 8,
        "b": 6,
        "a": 2,
    },
    212: {
        "e": 24,
        "d": 12,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    213: {
        "e": 24,
        "d": 12,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    214: {
        "i": 48,
        "h": 24,
        "g": 24,
        "f": 24,
        "e": 16,
        "d": 12,
        "c": 12,
        "b": 8,
        "a": 8,
    },
    215: {
        "j": 24,
        "i": 12,
        "h": 12,
        "g": 6,
        "f": 6,
        "e": 4,
        "d": 3,
        "c": 3,
        "b": 1,
        "a": 1,
    },
    216: {
        "i": 96,
        "h": 48,
        "g": 24,
        "f": 24,
        "e": 16,
        "d": 4,
        "c": 4,
        "b": 4,
        "a": 4,
    },
    217: {
        "h": 48,
        "g": 24,
        "f": 24,
        "e": 12,
        "d": 12,
        "c": 8,
        "b": 6,
        "a": 2,
    },
    218: {
        "i": 24,
        "h": 12,
        "g": 12,
        "f": 12,
        "e": 8,
        "d": 6,
        "c": 6,
        "b": 6,
        "a": 2,
    },
    219: {
        "h": 96,
        "g": 48,
        "f": 48,
        "e": 32,
        "d": 24,
        "c": 24,
        "b": 8,
        "a": 8,
    },
    220: {
        "e": 48,
        "d": 24,
        "c": 16,
        "b": 12,
        "a": 12,
    },
    221: {
        "n": 48,
        "m": 24,
        "l": 24,
        "k": 24,
        "j": 12,
        "i": 12,
        "h": 12,
        "g": 8,
        "f": 6,
        "e": 6,
        "d": 3,
        "c": 3,
        "b": 1,
        "a": 1,
    },
    222: {
        "i": 48,
        "h": 24,
        "g": 24,
        "f": 16,
        "e": 12,
        "d": 12,
        "c": 8,
        "b": 6,
        "a": 2,
    },
    223: {
        "l": 48,
        "k": 24,
        "j": 24,
        "i": 16,
        "h": 12,
        "g": 12,
        "f": 12,
        "e": 8,
        "d": 6,
        "c": 6,
        "b": 6,
        "a": 2,
    },
    224: {
        "l": 48,
        "k": 24,
        "j": 24,
        "i": 24,
        "h": 24,
        "g": 12,
        "f": 12,
        "e": 8,
        "d": 6,
        "c": 4,
        "b": 4,
        "a": 2,
    },
    225: {
        "l": 192,
        "k": 96,
        "j": 96,
        "i": 48,
        "h": 48,
        "g": 48,
        "f": 32,
        "e": 24,
        "d": 24,
        "c": 8,
        "b": 4,
        "a": 4,
    },
    226: {
        "j": 192,
        "i": 96,
        "h": 96,
        "g": 64,
        "f": 48,
        "e": 48,
        "d": 24,
        "c": 24,
        "b": 8,
        "a": 8,
    },
    227: {
        "i": 192,
        "h": 96,
        "g": 96,
        "f": 48,
        "e": 32,
        "d": 16,
        "c": 16,
        "b": 8,
        "a": 8,
    },
    228: {
        "h": 192,
        "g": 96,
        "f": 96,
        "e": 64,
        "d": 48,
        "c": 32,
        "b": 32,
        "a": 16,
    },
    229: {
        "l": 96,
        "k": 48,
        "j": 48,
        "i": 48,
        "h": 24,
        "g": 24,
        "f": 16,
        "e": 12,
        "d": 12,
        "c": 8,
        "b": 6,
        "a": 2,
    },
    230: {
        "h": 96,
        "g": 48,
        "f": 48,
        "e": 32,
        "d": 24,
        "c": 24,
        "b": 16,
        "a": 16,
    },
}

CENTERING_DIVISORS = {
     'P': 1,
     'C': 2,
     'I': 2,
     'F': 4,
     'R': 3,
}

# This is the "Transformed WP" column of the tables at the bottom of the page for each space group from https://cryst.ehu.es/cryst/get_set.html
WYCK_POS_XFORM_UNDER_NORMALIZER = {
    "1": [
        ["a"]
    ],
    "2": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "g", "f", "h", "d", "c", "e", "i"],
        ["c", "g", "a", "e", "d", "h", "b", "f", "i"],
        ["d", "f", "e", "a", "c", "b", "h", "g", "i"],
        ["e", "h", "d", "c", "a", "g", "f", "b", "i"],
        ["f", "d", "h", "b", "g", "a", "e", "c", "i"],
        ["g", "c", "b", "h", "f", "e", "a", "d", "i"],
        ["h", "e", "f", "g", "b", "c", "d", "a", "i"]
    ],
    "3": [
        ["a", "b", "c", "d", "e"],
        ["b", "a", "d", "c", "e"],
        ["c", "d", "a", "b", "e"],
        ["d", "c", "b", "a", "e"]
    ],
    "4": [
        ["a"]
    ],
    "5": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "6": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "7": [
        ["a"]
    ],
    "8": [
        ["a", "b"]
    ],
    "9": [
        ["a"]
    ],
    "10": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o"],
        ["b", "a", "f", "e", "d", "c", "h", "g", "i", "j", "k", "l", "n", "m", "o"],
        ["c", "f", "a", "g", "h", "b", "d", "e", "k", "l", "i", "j", "m", "n", "o"],
        ["d", "e", "g", "a", "b", "h", "c", "f", "j", "i", "l", "k", "m", "n", "o"],
        ["e", "d", "h", "b", "a", "g", "f", "c", "j", "i", "l", "k", "n", "m", "o"],
        ["f", "c", "b", "h", "g", "a", "e", "d", "k", "l", "i", "j", "n", "m", "o"],
        ["g", "h", "d", "c", "f", "e", "a", "b", "l", "k", "j", "i", "m", "n", "o"],
        ["h", "g", "e", "f", "c", "d", "b", "a", "l", "k", "j", "i", "n", "m", "o"]
    ],
    "11": [
        ["a", "b", "c", "d", "e", "f"],
        ["b", "a", "d", "c", "e", "f"],
        ["c", "d", "a", "b", "e", "f"],
        ["d", "c", "b", "a", "e", "f"]
    ],
    "12": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "i", "j"],
        ["c", "d", "a", "b", "f", "e", "h", "g", "i", "j"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i", "j"]
    ],
    "13": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["b", "a", "d", "c", "f", "e", "g"],
        ["c", "d", "a", "b", "e", "f", "g"],
        ["d", "c", "b", "a", "f", "e", "g"]
    ],
    "14": [
        ["a", "b", "c", "d", "e"],
        ["b", "a", "d", "c", "e"],
        ["c", "d", "a", "b", "e"],
        ["d", "c", "b", "a", "e"]
    ],
    "15": [
        ["a", "b", "c", "d", "e", "f"],
        ["a", "b", "d", "c", "e", "f"],
        ["b", "a", "c", "d", "e", "f"],
        ["b", "a", "d", "c", "e", "f"],
        ["c", "d", "a", "b", "e", "f"],
        ["c", "d", "b", "a", "e", "f"],
        ["d", "c", "a", "b", "e", "f"],
        ["d", "c", "b", "a", "e", "f"]
    ],
    "16": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u"],
        ["a", "b", "d", "c", "f", "e", "g", "h", "i", "k", "j", "l", "q", "s", "r", "t", "m", "o", "n", "p", "u"],
        ["a", "c", "b", "d", "e", "g", "f", "h", "m", "n", "o", "p", "i", "j", "k", "l", "q", "s", "r", "t", "u"],
        ["a", "c", "d", "b", "g", "e", "f", "h", "m", "o", "n", "p", "q", "r", "s", "t", "i", "k", "j", "l", "u"],
        ["a", "d", "b", "c", "f", "g", "e", "h", "q", "s", "r", "t", "i", "k", "j", "l", "m", "n", "o", "p", "u"],
        ["a", "d", "c", "b", "g", "f", "e", "h", "q", "r", "s", "t", "m", "o", "n", "p", "i", "j", "k", "l", "u"],
        ["b", "a", "e", "f", "c", "d", "h", "g", "i", "j", "k", "l", "o", "p", "m", "n", "r", "q", "t", "s", "u"],
        ["b", "a", "f", "e", "d", "c", "h", "g", "i", "k", "j", "l", "r", "t", "q", "s", "o", "m", "p", "n", "u"],
        ["b", "e", "a", "f", "c", "h", "d", "g", "o", "p", "m", "n", "i", "j", "k", "l", "r", "t", "q", "s", "u"],
        ["b", "e", "f", "a", "h", "c", "d", "g", "o", "m", "p", "n", "r", "q", "t", "s", "i", "k", "j", "l", "u"],
        ["b", "f", "a", "e", "d", "h", "c", "g", "r", "t", "q", "s", "i", "k", "j", "l", "o", "p", "m", "n", "u"],
        ["b", "f", "e", "a", "h", "d", "c", "g", "r", "q", "t", "s", "o", "m", "p", "n", "i", "j", "k", "l", "u"],
        ["c", "a", "e", "g", "b", "d", "h", "f", "m", "n", "o", "p", "k", "l", "i", "j", "s", "q", "t", "r", "u"],
        ["c", "a", "g", "e", "d", "b", "h", "f", "m", "o", "n", "p", "s", "t", "q", "r", "k", "i", "l", "j", "u"],
        ["c", "e", "a", "g", "b", "h", "d", "f", "k", "l", "i", "j", "m", "n", "o", "p", "s", "t", "q", "r", "u"],
        ["c", "e", "g", "a", "h", "b", "d", "f", "k", "i", "l", "j", "s", "q", "t", "r", "m", "o", "n", "p", "u"],
        ["c", "g", "a", "e", "d", "h", "b", "f", "s", "t", "q", "r", "m", "o", "n", "p", "k", "l", "i", "j", "u"],
        ["c", "g", "e", "a", "h", "d", "b", "f", "s", "q", "t", "r", "k", "i", "l", "j", "m", "n", "o", "p", "u"],
        ["d", "a", "f", "g", "b", "c", "h", "e", "q", "s", "r", "t", "j", "l", "i", "k", "n", "m", "p", "o", "u"],
        ["d", "a", "g", "f", "c", "b", "h", "e", "q", "r", "s", "t", "n", "p", "m", "o", "j", "i", "l", "k", "u"],
        ["d", "f", "a", "g", "b", "h", "c", "e", "j", "l", "i", "k", "q", "s", "r", "t", "n", "p", "m", "o", "u"],
        ["d", "f", "g", "a", "h", "b", "c", "e", "j", "i", "l", "k", "n", "m", "p", "o", "q", "r", "s", "t", "u"],
        ["d", "g", "a", "f", "c", "h", "b", "e", "n", "p", "m", "o", "q", "r", "s", "t", "j", "l", "i", "k", "u"],
        ["d", "g", "f", "a", "h", "c", "b", "e", "n", "m", "p", "o", "j", "i", "l", "k", "q", "s", "r", "t", "u"],
        ["e", "b", "c", "h", "a", "f", "g", "d", "o", "p", "m", "n", "k", "l", "i", "j", "t", "r", "s", "q", "u"],
        ["e", "b", "h", "c", "f", "a", "g", "d", "o", "m", "p", "n", "t", "s", "r", "q", "k", "i", "l", "j", "u"],
        ["e", "c", "b", "h", "a", "g", "f", "d", "k", "l", "i", "j", "o", "p", "m", "n", "t", "s", "r", "q", "u"],
        ["e", "c", "h", "b", "g", "a", "f", "d", "k", "i", "l", "j", "t", "r", "s", "q", "o", "m", "p", "n", "u"],
        ["e", "h", "b", "c", "f", "g", "a", "d", "t", "s", "r", "q", "o", "m", "p", "n", "k", "l", "i", "j", "u"],
        ["e", "h", "c", "b", "g", "f", "a", "d", "t", "r", "s", "q", "k", "i", "l", "j", "o", "p", "m", "n", "u"],
        ["f", "b", "d", "h", "a", "e", "g", "c", "r", "t", "q", "s", "j", "l", "i", "k", "p", "o", "n", "m", "u"],
        ["f", "b", "h", "d", "e", "a", "g", "c", "r", "q", "t", "s", "p", "n", "o", "m", "j", "i", "l", "k", "u"],
        ["f", "d", "b", "h", "a", "g", "e", "c", "j", "l", "i", "k", "r", "t", "q", "s", "p", "n", "o", "m", "u"],
        ["f", "d", "h", "b", "g", "a", "e", "c", "j", "i", "l", "k", "p", "o", "n", "m", "r", "q", "t", "s", "u"],
        ["f", "h", "b", "d", "e", "g", "a", "c", "p", "n", "o", "m", "r", "q", "t", "s", "j", "l", "i", "k", "u"],
        ["f", "h", "d", "b", "g", "e", "a", "c", "p", "o", "n", "m", "j", "i", "l", "k", "r", "t", "q", "s", "u"],
        ["g", "c", "d", "h", "a", "e", "f", "b", "s", "t", "q", "r", "n", "p", "m", "o", "l", "k", "j", "i", "u"],
        ["g", "c", "h", "d", "e", "a", "f", "b", "s", "q", "t", "r", "l", "j", "k", "i", "n", "m", "p", "o", "u"],
        ["g", "d", "c", "h", "a", "f", "e", "b", "n", "p", "m", "o", "s", "t", "q", "r", "l", "j", "k", "i", "u"],
        ["g", "d", "h", "c", "f", "a", "e", "b", "n", "m", "p", "o", "l", "k", "j", "i", "s", "q", "t", "r", "u"],
        ["g", "h", "c", "d", "e", "f", "a", "b", "l", "j", "k", "i", "s", "q", "t", "r", "n", "p", "m", "o", "u"],
        ["g", "h", "d", "c", "f", "e", "a", "b", "l", "k", "j", "i", "n", "m", "p", "o", "s", "t", "q", "r", "u"],
        ["h", "e", "f", "g", "b", "c", "d", "a", "t", "s", "r", "q", "p", "n", "o", "m", "l", "k", "j", "i", "u"],
        ["h", "e", "g", "f", "c", "b", "d", "a", "t", "r", "s", "q", "l", "j", "k", "i", "p", "o", "n", "m", "u"],
        ["h", "f", "e", "g", "b", "d", "c", "a", "p", "n", "o", "m", "t", "s", "r", "q", "l", "j", "k", "i", "u"],
        ["h", "f", "g", "e", "d", "b", "c", "a", "p", "o", "n", "m", "l", "k", "j", "i", "t", "r", "s", "q", "u"],
        ["h", "g", "e", "f", "c", "d", "b", "a", "l", "j", "k", "i", "t", "r", "s", "q", "p", "n", "o", "m", "u"],
        ["h", "g", "f", "e", "d", "c", "b", "a", "l", "k", "j", "i", "p", "o", "n", "m", "t", "s", "r", "q", "u"]
    ],
    "17": [
        ["a", "b", "c", "d", "e"],
        ["a", "b", "d", "c", "e"],
        ["b", "a", "c", "d", "e"],
        ["b", "a", "d", "c", "e"],
        ["c", "d", "a", "b", "e"],
        ["c", "d", "b", "a", "e"],
        ["d", "c", "a", "b", "e"],
        ["d", "c", "b", "a", "e"]
    ],
    "18": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "19": [
        ["a"]
    ],
    "20": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "21": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["a", "b", "c", "d", "g", "h", "e", "f", "i", "j", "k", "l"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "j", "i", "k", "l"],
        ["b", "a", "d", "c", "g", "h", "e", "f", "j", "i", "k", "l"],
        ["c", "d", "a", "b", "f", "e", "h", "g", "j", "i", "k", "l"],
        ["c", "d", "a", "b", "h", "g", "f", "e", "j", "i", "k", "l"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i", "j", "k", "l"],
        ["d", "c", "b", "a", "h", "g", "f", "e", "i", "j", "k", "l"]
    ],
    "22": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["a", "b", "c", "d", "e", "g", "f", "i", "h", "j", "k"],
        ["a", "b", "c", "d", "f", "e", "g", "h", "j", "i", "k"],
        ["a", "b", "c", "d", "f", "g", "e", "j", "h", "i", "k"],
        ["a", "b", "c", "d", "g", "e", "f", "i", "j", "h", "k"],
        ["a", "b", "c", "d", "g", "f", "e", "j", "i", "h", "k"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i", "j", "k"],
        ["a", "b", "d", "c", "e", "g", "f", "i", "h", "j", "k"],
        ["a", "b", "d", "c", "f", "e", "g", "h", "j", "i", "k"],
        ["a", "b", "d", "c", "f", "g", "e", "j", "h", "i", "k"],
        ["a", "b", "d", "c", "g", "e", "f", "i", "j", "h", "k"],
        ["a", "b", "d", "c", "g", "f", "e", "j", "i", "h", "k"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["b", "a", "c", "d", "e", "g", "f", "i", "h", "j", "k"],
        ["b", "a", "c", "d", "f", "e", "g", "h", "j", "i", "k"],
        ["b", "a", "c", "d", "f", "g", "e", "j", "h", "i", "k"],
        ["b", "a", "c", "d", "g", "e", "f", "i", "j", "h", "k"],
        ["b", "a", "c", "d", "g", "f", "e", "j", "i", "h", "k"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "i", "j", "k"],
        ["b", "a", "d", "c", "e", "g", "f", "i", "h", "j", "k"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "j", "i", "k"],
        ["b", "a", "d", "c", "f", "g", "e", "j", "h", "i", "k"],
        ["b", "a", "d", "c", "g", "e", "f", "i", "j", "h", "k"],
        ["b", "a", "d", "c", "g", "f", "e", "j", "i", "h", "k"],
        ["c", "d", "a", "b", "h", "i", "j", "e", "f", "g", "k"],
        ["c", "d", "a", "b", "h", "j", "i", "f", "e", "g", "k"],
        ["c", "d", "a", "b", "i", "h", "j", "e", "g", "f", "k"],
        ["c", "d", "a", "b", "i", "j", "h", "g", "e", "f", "k"],
        ["c", "d", "a", "b", "j", "h", "i", "f", "g", "e", "k"],
        ["c", "d", "a", "b", "j", "i", "h", "g", "f", "e", "k"],
        ["c", "d", "b", "a", "h", "i", "j", "e", "f", "g", "k"],
        ["c", "d", "b", "a", "h", "j", "i", "f", "e", "g", "k"],
        ["c", "d", "b", "a", "i", "h", "j", "e", "g", "f", "k"],
        ["c", "d", "b", "a", "i", "j", "h", "g", "e", "f", "k"],
        ["c", "d", "b", "a", "j", "h", "i", "f", "g", "e", "k"],
        ["c", "d", "b", "a", "j", "i", "h", "g", "f", "e", "k"],
        ["d", "c", "a", "b", "h", "i", "j", "e", "f", "g", "k"],
        ["d", "c", "a", "b", "h", "j", "i", "f", "e", "g", "k"],
        ["d", "c", "a", "b", "i", "h", "j", "e", "g", "f", "k"],
        ["d", "c", "a", "b", "i", "j", "h", "g", "e", "f", "k"],
        ["d", "c", "a", "b", "j", "h", "i", "f", "g", "e", "k"],
        ["d", "c", "a", "b", "j", "i", "h", "g", "f", "e", "k"],
        ["d", "c", "b", "a", "h", "i", "j", "e", "f", "g", "k"],
        ["d", "c", "b", "a", "h", "j", "i", "f", "e", "g", "k"],
        ["d", "c", "b", "a", "i", "h", "j", "e", "g", "f", "k"],
        ["d", "c", "b", "a", "i", "j", "h", "g", "e", "f", "k"],
        ["d", "c", "b", "a", "j", "h", "i", "f", "g", "e", "k"],
        ["d", "c", "b", "a", "j", "i", "h", "g", "f", "e", "k"]
    ],
    "23": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["a", "b", "d", "c", "e", "f", "i", "j", "g", "h", "k"],
        ["a", "c", "b", "d", "i", "j", "g", "h", "e", "f", "k"],
        ["a", "c", "d", "b", "i", "j", "e", "f", "g", "h", "k"],
        ["a", "d", "b", "c", "g", "h", "i", "j", "e", "f", "k"],
        ["a", "d", "c", "b", "g", "h", "e", "f", "i", "j", "k"],
        ["b", "a", "c", "d", "e", "f", "j", "i", "h", "g", "k"],
        ["b", "a", "d", "c", "e", "f", "h", "g", "j", "i", "k"],
        ["b", "c", "a", "d", "h", "g", "j", "i", "e", "f", "k"],
        ["b", "c", "d", "a", "h", "g", "e", "f", "j", "i", "k"],
        ["b", "d", "a", "c", "j", "i", "h", "g", "e", "f", "k"],
        ["b", "d", "c", "a", "j", "i", "e", "f", "h", "g", "k"],
        ["c", "a", "b", "d", "i", "j", "f", "e", "h", "g", "k"],
        ["c", "a", "d", "b", "i", "j", "h", "g", "f", "e", "k"],
        ["c", "b", "a", "d", "h", "g", "f", "e", "i", "j", "k"],
        ["c", "b", "d", "a", "h", "g", "i", "j", "f", "e", "k"],
        ["c", "d", "a", "b", "f", "e", "h", "g", "i", "j", "k"],
        ["c", "d", "b", "a", "f", "e", "i", "j", "h", "g", "k"],
        ["d", "a", "b", "c", "g", "h", "f", "e", "j", "i", "k"],
        ["d", "a", "c", "b", "g", "h", "j", "i", "f", "e", "k"],
        ["d", "b", "a", "c", "j", "i", "f", "e", "g", "h", "k"],
        ["d", "b", "c", "a", "j", "i", "g", "h", "f", "e", "k"],
        ["d", "c", "a", "b", "f", "e", "j", "i", "g", "h", "k"],
        ["d", "c", "b", "a", "f", "e", "g", "h", "j", "i", "k"]
    ],
    "24": [
        ["a", "b", "c", "d"],
        ["a", "c", "b", "d"],
        ["b", "a", "c", "d"],
        ["b", "c", "a", "d"],
        ["c", "a", "b", "d"],
        ["c", "b", "a", "d"]
    ],
    "25": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["a", "c", "b", "d", "g", "h", "e", "f", "i"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "i"],
        ["b", "d", "a", "c", "g", "h", "f", "e", "i"],
        ["c", "a", "d", "b", "h", "g", "e", "f", "i"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "i"],
        ["d", "b", "c", "a", "h", "g", "f", "e", "i"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i"]
    ],
    "26": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "27": [
        ["a", "b", "c", "d", "e"],
        ["a", "c", "b", "d", "e"],
        ["b", "a", "d", "c", "e"],
        ["b", "d", "a", "c", "e"],
        ["c", "a", "d", "b", "e"],
        ["c", "d", "a", "b", "e"],
        ["d", "b", "c", "a", "e"],
        ["d", "c", "b", "a", "e"]
    ],
    "28": [
        ["a", "b", "c", "d"],
        ["b", "a", "c", "d"]
    ],
    "29": [
        ["a"]
    ],
    "30": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "31": [
        ["a", "b"]
    ],
    "32": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "33": [
        ["a"]
    ],
    "34": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "35": [
        ["a", "b", "c", "d", "e", "f"],
        ["a", "b", "c", "e", "d", "f"],
        ["b", "a", "c", "d", "e", "f"],
        ["b", "a", "c", "e", "d", "f"]
    ],
    "36": [
        ["a", "b"]
    ],
    "37": [
        ["a", "b", "c", "d"],
        ["b", "a", "c", "d"]
    ],
    "38": [
        ["a", "b", "c", "d", "e", "f"],
        ["b", "a", "c", "e", "d", "f"]
    ],
    "39": [
        ["a", "b", "c", "d"],
        ["b", "a", "c", "d"]
    ],
    "40": [
        ["a", "b", "c"]
    ],
    "41": [
        ["a", "b"]
    ],
    "42": [
        ["a", "b", "c", "d", "e"],
        ["a", "b", "d", "c", "e"]
    ],
    "43": [
        ["a", "b"]
    ],
    "44": [
        ["a", "b", "c", "d", "e"],
        ["a", "b", "d", "c", "e"],
        ["b", "a", "c", "d", "e"],
        ["b", "a", "d", "c", "e"]
    ],
    "45": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "46": [
        ["a", "b", "c"]
    ],
    "47": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "A"],
        ["a", "b", "e", "f", "c", "d", "g", "h", "i", "k", "j", "l", "q", "r", "s", "t", "m", "n", "o", "p", "u", "v", "y", "z", "w", "x", "A"],
        ["a", "c", "b", "d", "e", "g", "f", "h", "q", "s", "r", "t", "m", "o", "n", "p", "i", "k", "j", "l", "y", "z", "w", "x", "u", "v", "A"],
        ["a", "c", "e", "g", "b", "d", "f", "h", "q", "r", "s", "t", "i", "k", "j", "l", "m", "o", "n", "p", "y", "z", "u", "v", "w", "x", "A"],
        ["a", "e", "b", "f", "c", "g", "d", "h", "m", "o", "n", "p", "q", "s", "r", "t", "i", "j", "k", "l", "w", "x", "y", "z", "u", "v", "A"],
        ["a", "e", "c", "g", "b", "f", "d", "h", "m", "n", "o", "p", "i", "j", "k", "l", "q", "s", "r", "t", "w", "x", "u", "v", "y", "z", "A"],
        ["b", "a", "d", "c", "f", "e", "h", "g", "i", "j", "k", "l", "o", "p", "m", "n", "s", "t", "q", "r", "v", "u", "w", "x", "y", "z", "A"],
        ["b", "a", "f", "e", "d", "c", "h", "g", "i", "k", "j", "l", "s", "t", "q", "r", "o", "p", "m", "n", "v", "u", "y", "z", "w", "x", "A"],
        ["b", "d", "a", "c", "f", "h", "e", "g", "s", "q", "t", "r", "o", "m", "p", "n", "i", "k", "j", "l", "y", "z", "w", "x", "v", "u", "A"],
        ["b", "d", "f", "h", "a", "c", "e", "g", "s", "t", "q", "r", "i", "k", "j", "l", "o", "m", "p", "n", "y", "z", "v", "u", "w", "x", "A"],
        ["b", "f", "a", "e", "d", "h", "c", "g", "o", "m", "p", "n", "s", "q", "t", "r", "i", "j", "k", "l", "w", "x", "y", "z", "v", "u", "A"],
        ["b", "f", "d", "h", "a", "e", "c", "g", "o", "p", "m", "n", "i", "j", "k", "l", "s", "q", "t", "r", "w", "x", "v", "u", "y", "z", "A"],
        ["c", "a", "d", "b", "g", "e", "h", "f", "q", "s", "r", "t", "n", "p", "m", "o", "j", "l", "i", "k", "z", "y", "w", "x", "u", "v", "A"],
        ["c", "a", "g", "e", "d", "b", "h", "f", "q", "r", "s", "t", "j", "l", "i", "k", "n", "p", "m", "o", "z", "y", "u", "v", "w", "x", "A"],
        ["c", "d", "a", "b", "g", "h", "e", "f", "j", "i", "l", "k", "n", "m", "p", "o", "q", "r", "s", "t", "u", "v", "w", "x", "z", "y", "A"],
        ["c", "d", "g", "h", "a", "b", "e", "f", "j", "l", "i", "k", "q", "r", "s", "t", "n", "m", "p", "o", "u", "v", "z", "y", "w", "x", "A"],
        ["c", "g", "a", "e", "d", "h", "b", "f", "n", "m", "p", "o", "j", "i", "l", "k", "q", "s", "r", "t", "w", "x", "u", "v", "z", "y", "A"],
        ["c", "g", "d", "h", "a", "e", "b", "f", "n", "p", "m", "o", "q", "s", "r", "t", "j", "i", "l", "k", "w", "x", "z", "y", "u", "v", "A"],
        ["d", "b", "c", "a", "h", "f", "g", "e", "s", "q", "t", "r", "p", "n", "o", "m", "j", "l", "i", "k", "z", "y", "w", "x", "v", "u", "A"],
        ["d", "b", "h", "f", "c", "a", "g", "e", "s", "t", "q", "r", "j", "l", "i", "k", "p", "n", "o", "m", "z", "y", "v", "u", "w", "x", "A"],
        ["d", "c", "b", "a", "h", "g", "f", "e", "j", "i", "l", "k", "p", "o", "n", "m", "s", "t", "q", "r", "v", "u", "w", "x", "z", "y", "A"],
        ["d", "c", "h", "g", "b", "a", "f", "e", "j", "l", "i", "k", "s", "t", "q", "r", "p", "o", "n", "m", "v", "u", "z", "y", "w", "x", "A"],
        ["d", "h", "b", "f", "c", "g", "a", "e", "p", "o", "n", "m", "j", "i", "l", "k", "s", "q", "t", "r", "w", "x", "v", "u", "z", "y", "A"],
        ["d", "h", "c", "g", "b", "f", "a", "e", "p", "n", "o", "m", "s", "q", "t", "r", "j", "i", "l", "k", "w", "x", "z", "y", "v", "u", "A"],
        ["e", "a", "f", "b", "g", "c", "h", "d", "m", "o", "n", "p", "r", "t", "q", "s", "k", "l", "i", "j", "x", "w", "y", "z", "u", "v", "A"],
        ["e", "a", "g", "c", "f", "b", "h", "d", "m", "n", "o", "p", "k", "l", "i", "j", "r", "t", "q", "s", "x", "w", "u", "v", "y", "z", "A"],
        ["e", "f", "a", "b", "g", "h", "c", "d", "k", "i", "l", "j", "r", "q", "t", "s", "m", "n", "o", "p", "u", "v", "y", "z", "x", "w", "A"],
        ["e", "f", "g", "h", "a", "b", "c", "d", "k", "l", "i", "j", "m", "n", "o", "p", "r", "q", "t", "s", "u", "v", "x", "w", "y", "z", "A"],
        ["e", "g", "a", "c", "f", "h", "b", "d", "r", "q", "t", "s", "k", "i", "l", "j", "m", "o", "n", "p", "y", "z", "u", "v", "x", "w", "A"],
        ["e", "g", "f", "h", "a", "c", "b", "d", "r", "t", "q", "s", "m", "o", "n", "p", "k", "i", "l", "j", "y", "z", "x", "w", "u", "v", "A"],
        ["f", "b", "e", "a", "h", "d", "g", "c", "o", "m", "p", "n", "t", "r", "s", "q", "k", "l", "i", "j", "x", "w", "y", "z", "v", "u", "A"],
        ["f", "b", "h", "d", "e", "a", "g", "c", "o", "p", "m", "n", "k", "l", "i", "j", "t", "r", "s", "q", "x", "w", "v", "u", "y", "z", "A"],
        ["f", "e", "b", "a", "h", "g", "d", "c", "k", "i", "l", "j", "t", "s", "r", "q", "o", "p", "m", "n", "v", "u", "y", "z", "x", "w", "A"],
        ["f", "e", "h", "g", "b", "a", "d", "c", "k", "l", "i", "j", "o", "p", "m", "n", "t", "s", "r", "q", "v", "u", "x", "w", "y", "z", "A"],
        ["f", "h", "b", "d", "e", "g", "a", "c", "t", "s", "r", "q", "k", "i", "l", "j", "o", "m", "p", "n", "y", "z", "v", "u", "x", "w", "A"],
        ["f", "h", "e", "g", "b", "d", "a", "c", "t", "r", "s", "q", "o", "m", "p", "n", "k", "i", "l", "j", "y", "z", "x", "w", "v", "u", "A"],
        ["g", "c", "e", "a", "h", "d", "f", "b", "n", "m", "p", "o", "l", "k", "j", "i", "r", "t", "q", "s", "x", "w", "u", "v", "z", "y", "A"],
        ["g", "c", "h", "d", "e", "a", "f", "b", "n", "p", "m", "o", "r", "t", "q", "s", "l", "k", "j", "i", "x", "w", "z", "y", "u", "v", "A"],
        ["g", "e", "c", "a", "h", "f", "d", "b", "r", "q", "t", "s", "l", "j", "k", "i", "n", "p", "m", "o", "z", "y", "u", "v", "x", "w", "A"],
        ["g", "e", "h", "f", "c", "a", "d", "b", "r", "t", "q", "s", "n", "p", "m", "o", "l", "j", "k", "i", "z", "y", "x", "w", "u", "v", "A"],
        ["g", "h", "c", "d", "e", "f", "a", "b", "l", "j", "k", "i", "r", "q", "t", "s", "n", "m", "p", "o", "u", "v", "z", "y", "x", "w", "A"],
        ["g", "h", "e", "f", "c", "d", "a", "b", "l", "k", "j", "i", "n", "m", "p", "o", "r", "q", "t", "s", "u", "v", "x", "w", "z", "y", "A"],
        ["h", "d", "f", "b", "g", "c", "e", "a", "p", "o", "n", "m", "l", "k", "j", "i", "t", "r", "s", "q", "x", "w", "v", "u", "z", "y", "A"],
        ["h", "d", "g", "c", "f", "b", "e", "a", "p", "n", "o", "m", "t", "r", "s", "q", "l", "k", "j", "i", "x", "w", "z", "y", "v", "u", "A"],
        ["h", "f", "d", "b", "g", "e", "c", "a", "t", "s", "r", "q", "l", "j", "k", "i", "p", "n", "o", "m", "z", "y", "v", "u", "x", "w", "A"],
        ["h", "f", "g", "e", "d", "b", "c", "a", "t", "r", "s", "q", "p", "n", "o", "m", "l", "j", "k", "i", "z", "y", "x", "w", "v", "u", "A"],
        ["h", "g", "d", "c", "f", "e", "b", "a", "l", "j", "k", "i", "t", "s", "r", "q", "p", "o", "n", "m", "v", "u", "z", "y", "x", "w", "A"],
        ["h", "g", "f", "e", "d", "c", "b", "a", "l", "k", "j", "i", "p", "o", "n", "m", "t", "s", "r", "q", "v", "u", "x", "w", "z", "y", "A"]
    ],
    "48": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"],
        ["a", "b", "c", "d", "f", "e", "g", "h", "i", "j", "k", "l", "m"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "k", "l", "i", "j", "m"],
        ["a", "b", "d", "c", "f", "e", "g", "h", "k", "l", "i", "j", "m"],
        ["a", "c", "b", "d", "e", "f", "k", "l", "i", "j", "g", "h", "m"],
        ["a", "c", "b", "d", "f", "e", "k", "l", "i", "j", "g", "h", "m"],
        ["a", "c", "d", "b", "e", "f", "k", "l", "g", "h", "i", "j", "m"],
        ["a", "c", "d", "b", "f", "e", "k", "l", "g", "h", "i", "j", "m"],
        ["a", "d", "b", "c", "e", "f", "i", "j", "k", "l", "g", "h", "m"],
        ["a", "d", "b", "c", "f", "e", "i", "j", "k", "l", "g", "h", "m"],
        ["a", "d", "c", "b", "e", "f", "i", "j", "g", "h", "k", "l", "m"],
        ["a", "d", "c", "b", "f", "e", "i", "j", "g", "h", "k", "l", "m"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "l", "k", "j", "i", "m"],
        ["b", "a", "c", "d", "f", "e", "g", "h", "l", "k", "j", "i", "m"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "j", "i", "l", "k", "m"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "j", "i", "l", "k", "m"],
        ["b", "c", "a", "d", "e", "f", "j", "i", "l", "k", "g", "h", "m"],
        ["b", "c", "a", "d", "f", "e", "j", "i", "l", "k", "g", "h", "m"],
        ["b", "c", "d", "a", "e", "f", "j", "i", "g", "h", "l", "k", "m"],
        ["b", "c", "d", "a", "f", "e", "j", "i", "g", "h", "l", "k", "m"],
        ["b", "d", "a", "c", "e", "f", "l", "k", "j", "i", "g", "h", "m"],
        ["b", "d", "a", "c", "f", "e", "l", "k", "j", "i", "g", "h", "m"],
        ["b", "d", "c", "a", "e", "f", "l", "k", "g", "h", "j", "i", "m"],
        ["b", "d", "c", "a", "f", "e", "l", "k", "g", "h", "j", "i", "m"],
        ["c", "a", "b", "d", "e", "f", "k", "l", "h", "g", "j", "i", "m"],
        ["c", "a", "b", "d", "f", "e", "k", "l", "h", "g", "j", "i", "m"],
        ["c", "a", "d", "b", "e", "f", "k", "l", "j", "i", "h", "g", "m"],
        ["c", "a", "d", "b", "f", "e", "k", "l", "j", "i", "h", "g", "m"],
        ["c", "b", "a", "d", "e", "f", "j", "i", "h", "g", "k", "l", "m"],
        ["c", "b", "a", "d", "f", "e", "j", "i", "h", "g", "k", "l", "m"],
        ["c", "b", "d", "a", "e", "f", "j", "i", "k", "l", "h", "g", "m"],
        ["c", "b", "d", "a", "f", "e", "j", "i", "k", "l", "h", "g", "m"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "j", "i", "k", "l", "m"],
        ["c", "d", "a", "b", "f", "e", "h", "g", "j", "i", "k", "l", "m"],
        ["c", "d", "b", "a", "e", "f", "h", "g", "k", "l", "j", "i", "m"],
        ["c", "d", "b", "a", "f", "e", "h", "g", "k", "l", "j", "i", "m"],
        ["d", "a", "b", "c", "e", "f", "i", "j", "h", "g", "l", "k", "m"],
        ["d", "a", "b", "c", "f", "e", "i", "j", "h", "g", "l", "k", "m"],
        ["d", "a", "c", "b", "e", "f", "i", "j", "l", "k", "h", "g", "m"],
        ["d", "a", "c", "b", "f", "e", "i", "j", "l", "k", "h", "g", "m"],
        ["d", "b", "a", "c", "e", "f", "l", "k", "h", "g", "i", "j", "m"],
        ["d", "b", "a", "c", "f", "e", "l", "k", "h", "g", "i", "j", "m"],
        ["d", "b", "c", "a", "e", "f", "l", "k", "i", "j", "h", "g", "m"],
        ["d", "b", "c", "a", "f", "e", "l", "k", "i", "j", "h", "g", "m"],
        ["d", "c", "a", "b", "e", "f", "h", "g", "l", "k", "i", "j", "m"],
        ["d", "c", "a", "b", "f", "e", "h", "g", "l", "k", "i", "j", "m"],
        ["d", "c", "b", "a", "e", "f", "h", "g", "i", "j", "l", "k", "m"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i", "j", "l", "k", "m"]
    ],
    "49": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r"],
        ["a", "b", "d", "c", "e", "g", "f", "h", "k", "l", "i", "j", "m", "n", "p", "o", "q", "r"],
        ["b", "a", "c", "d", "h", "f", "g", "e", "l", "k", "j", "i", "n", "m", "o", "p", "q", "r"],
        ["b", "a", "d", "c", "h", "g", "f", "e", "j", "i", "l", "k", "n", "m", "p", "o", "q", "r"],
        ["c", "d", "a", "b", "g", "h", "e", "f", "j", "i", "k", "l", "o", "p", "m", "n", "q", "r"],
        ["c", "d", "b", "a", "g", "e", "h", "f", "k", "l", "j", "i", "o", "p", "n", "m", "q", "r"],
        ["d", "c", "a", "b", "f", "h", "e", "g", "l", "k", "i", "j", "p", "o", "m", "n", "q", "r"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i", "j", "l", "k", "p", "o", "n", "m", "q", "r"]
    ],
    "50": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"],
        ["a", "b", "c", "d", "e", "f", "i", "j", "g", "h", "k", "l", "m"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "i", "j", "l", "k", "m"],
        ["b", "a", "d", "c", "e", "f", "i", "j", "g", "h", "l", "k", "m"],
        ["c", "d", "a", "b", "f", "e", "h", "g", "j", "i", "l", "k", "m"],
        ["c", "d", "a", "b", "f", "e", "j", "i", "h", "g", "l", "k", "m"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "j", "i", "k", "l", "m"],
        ["d", "c", "b", "a", "f", "e", "j", "i", "h", "g", "k", "l", "m"]
    ],
    "51": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "j", "i", "k", "l"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "i", "j", "k", "l"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "j", "i", "k", "l"]
    ],
    "52": [
        ["a", "b", "c", "d", "e"],
        ["b", "a", "c", "d", "e"]
    ],
    "53": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "i"],
        ["c", "d", "a", "b", "f", "e", "g", "h", "i"],
        ["d", "c", "b", "a", "f", "e", "g", "h", "i"]
    ],
    "54": [
        ["a", "b", "c", "d", "e", "f"],
        ["b", "a", "c", "e", "d", "f"]
    ],
    "55": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "d", "c", "e", "f", "h", "g", "i"],
        ["c", "d", "a", "b", "f", "e", "g", "h", "i"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i"]
    ],
    "56": [
        ["a", "b", "c", "d", "e"],
        ["a", "b", "d", "c", "e"],
        ["b", "a", "c", "d", "e"],
        ["b", "a", "d", "c", "e"]
    ],
    "57": [
        ["a", "b", "c", "d", "e"],
        ["b", "a", "c", "d", "e"]
    ],
    "58": [
        ["a", "b", "c", "d", "e", "f", "g", "h"],
        ["a", "b", "d", "c", "e", "f", "g", "h"],
        ["b", "a", "c", "d", "e", "f", "g", "h"],
        ["b", "a", "d", "c", "e", "f", "g", "h"],
        ["c", "d", "a", "b", "f", "e", "g", "h"],
        ["c", "d", "b", "a", "f", "e", "g", "h"],
        ["d", "c", "a", "b", "f", "e", "g", "h"],
        ["d", "c", "b", "a", "f", "e", "g", "h"]
    ],
    "59": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["a", "b", "c", "d", "f", "e", "g"],
        ["a", "b", "d", "c", "e", "f", "g"],
        ["a", "b", "d", "c", "f", "e", "g"],
        ["b", "a", "c", "d", "e", "f", "g"],
        ["b", "a", "c", "d", "f", "e", "g"],
        ["b", "a", "d", "c", "e", "f", "g"],
        ["b", "a", "d", "c", "f", "e", "g"]
    ],
    "60": [
        ["a", "b", "c", "d"],
        ["b", "a", "c", "d"]
    ],
    "61": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "62": [
        ["a", "b", "c", "d"],
        ["b", "a", "c", "d"]
    ],
    "63": [
        ["a", "b", "c", "d", "e", "f", "g", "h"],
        ["b", "a", "c", "d", "e", "f", "g", "h"]
    ],
    "64": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["b", "a", "c", "d", "e", "f", "g"]
    ],
    "65": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r"],
        ["a", "b", "c", "d", "e", "f", "i", "j", "g", "h", "k", "l", "m", "o", "n", "p", "q", "r"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "i", "j", "l", "k", "m", "n", "o", "p", "q", "r"],
        ["b", "a", "d", "c", "e", "f", "i", "j", "g", "h", "l", "k", "m", "o", "n", "p", "q", "r"],
        ["c", "d", "a", "b", "f", "e", "h", "g", "j", "i", "l", "k", "m", "n", "o", "q", "p", "r"],
        ["c", "d", "a", "b", "f", "e", "j", "i", "h", "g", "l", "k", "m", "o", "n", "q", "p", "r"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "j", "i", "k", "l", "m", "n", "o", "q", "p", "r"],
        ["d", "c", "b", "a", "f", "e", "j", "i", "h", "g", "k", "l", "m", "o", "n", "q", "p", "r"]
    ],
    "66": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"],
        ["a", "b", "c", "d", "e", "f", "h", "g", "i", "j", "k", "l", "m"],
        ["a", "b", "c", "d", "f", "e", "g", "h", "i", "j", "k", "l", "m"],
        ["a", "b", "c", "d", "f", "e", "h", "g", "i", "j", "k", "l", "m"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "j", "i", "k", "l", "m"],
        ["b", "a", "d", "c", "e", "f", "h", "g", "j", "i", "k", "l", "m"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "j", "i", "k", "l", "m"],
        ["b", "a", "d", "c", "f", "e", "h", "g", "j", "i", "k", "l", "m"]
    ],
    "67": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o"],
        ["a", "b", "e", "f", "c", "d", "g", "j", "k", "h", "i", "l", "n", "m", "o"],
        ["b", "a", "d", "c", "f", "e", "g", "i", "h", "k", "j", "l", "m", "n", "o"],
        ["b", "a", "f", "e", "d", "c", "g", "k", "j", "i", "h", "l", "n", "m", "o"]
    ],
    "68": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["a", "b", "d", "c", "f", "e", "g", "h", "i"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "i"]
    ],
    "69": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p"],
        ["a", "b", "c", "e", "d", "f", "g", "i", "h", "k", "j", "l", "m", "o", "n", "p"],
        ["a", "b", "d", "c", "e", "f", "h", "g", "i", "j", "l", "k", "n", "m", "o", "p"],
        ["a", "b", "d", "e", "c", "f", "h", "i", "g", "l", "j", "k", "n", "o", "m", "p"],
        ["a", "b", "e", "c", "d", "f", "i", "g", "h", "k", "l", "j", "o", "m", "n", "p"],
        ["a", "b", "e", "d", "c", "f", "i", "h", "g", "l", "k", "j", "o", "n", "m", "p"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p"],
        ["b", "a", "c", "e", "d", "f", "g", "i", "h", "k", "j", "l", "m", "o", "n", "p"],
        ["b", "a", "d", "c", "e", "f", "h", "g", "i", "j", "l", "k", "n", "m", "o", "p"],
        ["b", "a", "d", "e", "c", "f", "h", "i", "g", "l", "j", "k", "n", "o", "m", "p"],
        ["b", "a", "e", "c", "d", "f", "i", "g", "h", "k", "l", "j", "o", "m", "n", "p"],
        ["b", "a", "e", "d", "c", "f", "i", "h", "g", "l", "k", "j", "o", "n", "m", "p"]
    ],
    "70": [
        ["a", "b", "c", "d", "e", "f", "g", "h"],
        ["a", "b", "c", "d", "e", "g", "f", "h"],
        ["a", "b", "c", "d", "f", "e", "g", "h"],
        ["a", "b", "c", "d", "f", "g", "e", "h"],
        ["a", "b", "c", "d", "g", "e", "f", "h"],
        ["a", "b", "c", "d", "g", "f", "e", "h"],
        ["b", "a", "d", "c", "e", "f", "g", "h"],
        ["b", "a", "d", "c", "e", "g", "f", "h"],
        ["b", "a", "d", "c", "f", "e", "g", "h"],
        ["b", "a", "d", "c", "f", "g", "e", "h"],
        ["b", "a", "d", "c", "g", "e", "f", "h"],
        ["b", "a", "d", "c", "g", "f", "e", "h"]
    ],
    "71": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o"],
        ["a", "b", "d", "c", "e", "f", "i", "j", "g", "h", "k", "l", "n", "m", "o"],
        ["a", "c", "b", "d", "i", "j", "g", "h", "e", "f", "k", "n", "m", "l", "o"],
        ["a", "c", "d", "b", "i", "j", "e", "f", "g", "h", "k", "n", "l", "m", "o"],
        ["a", "d", "b", "c", "g", "h", "i", "j", "e", "f", "k", "m", "n", "l", "o"],
        ["a", "d", "c", "b", "g", "h", "e", "f", "i", "j", "k", "m", "l", "n", "o"],
        ["b", "a", "c", "d", "e", "f", "j", "i", "h", "g", "k", "l", "n", "m", "o"],
        ["b", "a", "d", "c", "e", "f", "h", "g", "j", "i", "k", "l", "m", "n", "o"],
        ["b", "c", "a", "d", "h", "g", "j", "i", "e", "f", "k", "m", "n", "l", "o"],
        ["b", "c", "d", "a", "h", "g", "e", "f", "j", "i", "k", "m", "l", "n", "o"],
        ["b", "d", "a", "c", "j", "i", "h", "g", "e", "f", "k", "n", "m", "l", "o"],
        ["b", "d", "c", "a", "j", "i", "e", "f", "h", "g", "k", "n", "l", "m", "o"],
        ["c", "a", "b", "d", "i", "j", "f", "e", "h", "g", "k", "n", "l", "m", "o"],
        ["c", "a", "d", "b", "i", "j", "h", "g", "f", "e", "k", "n", "m", "l", "o"],
        ["c", "b", "a", "d", "h", "g", "f", "e", "i", "j", "k", "m", "l", "n", "o"],
        ["c", "b", "d", "a", "h", "g", "i", "j", "f", "e", "k", "m", "n", "l", "o"],
        ["c", "d", "a", "b", "f", "e", "h", "g", "i", "j", "k", "l", "m", "n", "o"],
        ["c", "d", "b", "a", "f", "e", "i", "j", "h", "g", "k", "l", "n", "m", "o"],
        ["d", "a", "b", "c", "g", "h", "f", "e", "j", "i", "k", "m", "l", "n", "o"],
        ["d", "a", "c", "b", "g", "h", "j", "i", "f", "e", "k", "m", "n", "l", "o"],
        ["d", "b", "a", "c", "j", "i", "f", "e", "g", "h", "k", "n", "l", "m", "o"],
        ["d", "b", "c", "a", "j", "i", "g", "h", "f", "e", "k", "n", "m", "l", "o"],
        ["d", "c", "a", "b", "f", "e", "j", "i", "g", "h", "k", "l", "n", "m", "o"],
        ["d", "c", "b", "a", "f", "e", "g", "h", "j", "i", "k", "l", "m", "n", "o"]
    ],
    "72": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["a", "b", "c", "d", "e", "g", "f", "h", "i", "j", "k"],
        ["b", "a", "d", "c", "e", "f", "g", "i", "h", "j", "k"],
        ["b", "a", "d", "c", "e", "g", "f", "i", "h", "j", "k"]
    ],
    "73": [
        ["a", "b", "c", "d", "e", "f"],
        ["a", "b", "d", "e", "c", "f"],
        ["a", "b", "e", "c", "d", "f"],
        ["b", "a", "c", "e", "d", "f"],
        ["b", "a", "d", "c", "e", "f"],
        ["b", "a", "e", "d", "c", "f"]
    ],
    "74": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i", "j"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "i", "j"],
        ["c", "d", "a", "b", "e", "g", "f", "i", "h", "j"],
        ["c", "d", "b", "a", "e", "g", "f", "i", "h", "j"],
        ["d", "c", "a", "b", "e", "g", "f", "i", "h", "j"],
        ["d", "c", "b", "a", "e", "g", "f", "i", "h", "j"]
    ],
    "75": [
        ["a", "b", "c", "d"],
        ["b", "a", "c", "d"]
    ],
    "76": [
        ["a"]
    ],
    "77": [
        ["a", "b", "c", "d"],
        ["b", "a", "c", "d"]
    ],
    "78": [
        ["a"]
    ],
    "79": [
        ["a", "b", "c"]
    ],
    "80": [
        ["a", "b"]
    ],
    "81": [
        ["a", "b", "c", "d", "e", "f", "g", "h"],
        ["b", "a", "d", "c", "e", "f", "g", "h"],
        ["c", "d", "a", "b", "f", "e", "g", "h"],
        ["d", "c", "b", "a", "f", "e", "g", "h"]
    ],
    "82": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["a", "b", "d", "c", "e", "f", "g"],
        ["b", "a", "c", "d", "e", "f", "g"],
        ["b", "a", "d", "c", "e", "f", "g"],
        ["c", "d", "a", "b", "f", "e", "g"],
        ["c", "d", "b", "a", "f", "e", "g"],
        ["d", "c", "a", "b", "f", "e", "g"],
        ["d", "c", "b", "a", "f", "e", "g"]
    ],
    "83": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "i", "k", "j", "l"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "i", "j", "k", "l"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i", "k", "j", "l"]
    ],
    "84": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i", "j", "k"],
        ["b", "a", "c", "d", "f", "e", "h", "g", "i", "j", "k"],
        ["b", "a", "d", "c", "f", "e", "h", "g", "i", "j", "k"]
    ],
    "85": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["b", "a", "c", "e", "d", "f", "g"]
    ],
    "86": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["a", "b", "d", "c", "e", "f", "g"],
        ["b", "a", "c", "d", "e", "f", "g"],
        ["b", "a", "d", "c", "e", "f", "g"]
    ],
    "87": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "i"]
    ],
    "88": [
        ["a", "b", "c", "d", "e", "f"],
        ["b", "a", "d", "c", "e", "f"]
    ],
    "89": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "i", "k", "j", "n", "o", "l", "m", "p"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "i", "j", "k", "o", "n", "m", "l", "p"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i", "k", "j", "m", "l", "o", "n", "p"]
    ],
    "90": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["b", "a", "c", "d", "f", "e", "g"]
    ],
    "91": [
        ["a", "b", "c", "d"],
        ["b", "a", "c", "d"]
    ],
    "92": [
        ["a", "b"]
    ],
    "93": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p"],
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "o", "n", "p"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i", "l", "m", "j", "k", "n", "o", "p"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i", "l", "m", "j", "k", "o", "n", "p"],
        ["b", "a", "c", "d", "f", "e", "h", "g", "i", "k", "j", "m", "l", "n", "o", "p"],
        ["b", "a", "c", "d", "f", "e", "h", "g", "i", "k", "j", "m", "l", "o", "n", "p"],
        ["b", "a", "d", "c", "f", "e", "h", "g", "i", "m", "l", "k", "j", "n", "o", "p"],
        ["b", "a", "d", "c", "f", "e", "h", "g", "i", "m", "l", "k", "j", "o", "n", "p"]
    ],
    "94": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["a", "b", "c", "d", "f", "e", "g"],
        ["b", "a", "c", "d", "e", "f", "g"],
        ["b", "a", "c", "d", "f", "e", "g"]
    ],
    "95": [
        ["a", "b", "c", "d"],
        ["b", "a", "c", "d"]
    ],
    "96": [
        ["a", "b"]
    ],
    "97": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["b", "a", "c", "d", "e", "f", "g", "i", "h", "j", "k"]
    ],
    "98": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["a", "b", "c", "e", "d", "f", "g"],
        ["b", "a", "c", "d", "e", "f", "g"],
        ["b", "a", "c", "e", "d", "f", "g"]
    ],
    "99": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["b", "a", "c", "d", "f", "e", "g"]
    ],
    "100": [
        ["a", "b", "c", "d"]
    ],
    "101": [
        ["a", "b", "c", "d", "e"],
        ["b", "a", "c", "d", "e"]
    ],
    "102": [
        ["a", "b", "c", "d"]
    ],
    "103": [
        ["a", "b", "c", "d"],
        ["b", "a", "c", "d"]
    ],
    "104": [
        ["a", "b", "c"]
    ],
    "105": [
        ["a", "b", "c", "d", "e", "f"],
        ["b", "a", "c", "e", "d", "f"]
    ],
    "106": [
        ["a", "b", "c"]
    ],
    "107": [
        ["a", "b", "c", "d", "e"]
    ],
    "108": [
        ["a", "b", "c", "d"]
    ],
    "109": [
        ["a", "b", "c"]
    ],
    "110": [
        ["a", "b"]
    ],
    "111": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o"],
        ["b", "a", "d", "c", "f", "e", "h", "g", "j", "i", "l", "k", "m", "n", "o"],
        ["c", "d", "a", "b", "f", "e", "g", "h", "k", "l", "i", "j", "m", "n", "o"],
        ["d", "c", "b", "a", "e", "f", "h", "g", "l", "k", "j", "i", "m", "n", "o"]
    ],
    "112": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"],
        ["a", "d", "c", "b", "e", "f", "j", "i", "h", "g", "k", "l", "m", "n"],
        ["c", "b", "a", "d", "f", "e", "h", "g", "j", "i", "l", "k", "m", "n"],
        ["c", "d", "a", "b", "f", "e", "i", "j", "g", "h", "l", "k", "m", "n"]
    ],
    "113": [
        ["a", "b", "c", "d", "e", "f"],
        ["b", "a", "c", "d", "e", "f"]
    ],
    "114": [
        ["a", "b", "c", "d", "e"],
        ["b", "a", "c", "d", "e"]
    ],
    "115": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "i", "k", "j", "l"],
        ["c", "d", "a", "b", "f", "e", "g", "i", "h", "k", "j", "l"],
        ["d", "c", "b", "a", "e", "f", "g", "i", "h", "j", "k", "l"]
    ],
    "116": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["a", "b", "c", "d", "f", "e", "g", "h", "i", "j"],
        ["b", "a", "d", "c", "e", "f", "h", "g", "i", "j"],
        ["b", "a", "d", "c", "f", "e", "h", "g", "i", "j"]
    ],
    "117": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "d", "c", "e", "f", "h", "g", "i"]
    ],
    "118": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["a", "b", "c", "d", "e", "g", "f", "h", "i"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i"],
        ["a", "b", "d", "c", "e", "g", "f", "h", "i"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "c", "d", "e", "g", "f", "h", "i"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "i"],
        ["b", "a", "d", "c", "e", "g", "f", "h", "i"]
    ],
    "119": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i", "j"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "i", "j"],
        ["c", "d", "a", "b", "f", "e", "h", "g", "i", "j"],
        ["c", "d", "b", "a", "f", "e", "h", "g", "i", "j"],
        ["d", "c", "a", "b", "f", "e", "h", "g", "i", "j"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i", "j"]
    ],
    "120": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["d", "c", "b", "a", "h", "g", "f", "e", "i"]
    ],
    "121": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["b", "a", "c", "d", "e", "g", "f", "h", "i", "j"]
    ],
    "122": [
        ["a", "b", "c", "d", "e"],
        ["b", "a", "c", "d", "e"]
    ],
    "123": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "i", "k", "j", "m", "l", "o", "n", "q", "p", "r", "s", "t", "u"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "i", "j", "k", "n", "o", "l", "m", "p", "q", "r", "t", "s", "u"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i", "k", "j", "o", "n", "m", "l", "q", "p", "r", "t", "s", "u"]
    ],
    "124": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "i", "j", "l", "k", "m", "n"]
    ],
    "125": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "j", "i", "l", "k", "m", "n"]
    ],
    "126": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "j", "i", "k"]
    ],
    "127": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["b", "a", "d", "c", "e", "f", "h", "g", "j", "i", "k", "l"]
    ],
    "128": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "i"]
    ],
    "129": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["b", "a", "c", "e", "d", "f", "h", "g", "i", "j", "k"]
    ],
    "130": [
        ["a", "b", "c", "d", "e", "f", "g"]
    ],
    "131": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i", "l", "m", "j", "k", "n", "o", "p", "q", "r"],
        ["b", "a", "c", "d", "f", "e", "h", "g", "i", "k", "j", "m", "l", "n", "p", "o", "q", "r"],
        ["b", "a", "d", "c", "f", "e", "h", "g", "i", "m", "l", "k", "j", "n", "p", "o", "q", "r"]
    ],
    "132": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p"],
        ["a", "b", "c", "d", "e", "f", "g", "h", "j", "i", "k", "l", "m", "n", "o", "p"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "i", "j", "k", "m", "l", "n", "o", "p"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "j", "i", "k", "m", "l", "n", "o", "p"]
    ],
    "133": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["a", "b", "c", "d", "e", "f", "g", "i", "h", "j", "k"]
    ],
    "134": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"],
        ["a", "b", "c", "d", "f", "e", "g", "h", "i", "j", "l", "k", "m", "n"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "j", "i", "k", "l", "m", "n"],
        ["b", "a", "c", "d", "f", "e", "g", "h", "j", "i", "l", "k", "m", "n"]
    ],
    "135": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
    ],
    "136": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["a", "b", "c", "d", "e", "g", "f", "h", "i", "j", "k"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["b", "a", "c", "d", "e", "g", "f", "h", "i", "j", "k"]
    ],
    "137": [
        ["a", "b", "c", "d", "e", "f", "g", "h"],
        ["b", "a", "c", "d", "e", "f", "g", "h"]
    ],
    "138": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["a", "b", "d", "c", "e", "f", "h", "g", "i", "j"]
    ],
    "139": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "j", "i", "k", "l", "m", "n", "o"]
    ],
    "140": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    ],
    "141": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "i"]
    ],
    "142": [
        ["a", "b", "c", "d", "e", "f", "g"]
    ],
    "143": [
        ["a", "b", "c", "d"],
        ["a", "c", "b", "d"],
        ["b", "a", "c", "d"],
        ["b", "c", "a", "d"],
        ["c", "a", "b", "d"],
        ["c", "b", "a", "d"]
    ],
    "144": [
        ["a"]
    ],
    "145": [
        ["a"]
    ],
    "146": [
        ["a", "b"]
    ],
    "147": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["b", "a", "c", "d", "f", "e", "g"]
    ],
    "148": [
        ["a", "b", "c", "d", "e", "f"],
        ["b", "a", "c", "e", "d", "f"]
    ],
    "149": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["a", "b", "e", "f", "c", "d", "g", "i", "h", "j", "k", "l"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "i", "k", "j", "l"],
        ["b", "a", "f", "e", "d", "c", "g", "i", "h", "k", "j", "l"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "i", "j", "k", "l"],
        ["c", "d", "e", "f", "a", "b", "h", "i", "g", "j", "k", "l"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i", "k", "j", "l"],
        ["d", "c", "f", "e", "b", "a", "h", "i", "g", "k", "j", "l"],
        ["e", "f", "a", "b", "c", "d", "i", "g", "h", "j", "k", "l"],
        ["e", "f", "c", "d", "a", "b", "i", "h", "g", "j", "k", "l"],
        ["f", "e", "b", "a", "d", "c", "i", "g", "h", "k", "j", "l"],
        ["f", "e", "d", "c", "b", "a", "i", "h", "g", "k", "j", "l"]
    ],
    "150": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["b", "a", "c", "d", "f", "e", "g"]
    ],
    "151": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "152": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "153": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "154": [
        ["a", "b", "c"],
        ["b", "a", "c"]
    ],
    "155": [
        ["a", "b", "c", "d", "e", "f"],
        ["b", "a", "c", "e", "d", "f"]
    ],
    "156": [
        ["a", "b", "c", "d", "e"],
        ["a", "c", "b", "d", "e"],
        ["b", "a", "c", "d", "e"],
        ["b", "c", "a", "d", "e"],
        ["c", "a", "b", "d", "e"],
        ["c", "b", "a", "d", "e"]
    ],
    "157": [
        ["a", "b", "c", "d"]
    ],
    "158": [
        ["a", "b", "c", "d"],
        ["a", "c", "b", "d"],
        ["b", "a", "c", "d"],
        ["b", "c", "a", "d"],
        ["c", "a", "b", "d"],
        ["c", "b", "a", "d"]
    ],
    "159": [
        ["a", "b", "c"]
    ],
    "160": [
        ["a", "b", "c"]
    ],
    "161": [
        ["a", "b"]
    ],
    "162": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["b", "a", "d", "c", "e", "g", "f", "h", "j", "i", "k", "l"]
    ],
    "163": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i"]
    ],
    "164": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["b", "a", "c", "d", "f", "e", "h", "g", "i", "j"]
    ],
    "165": [
        ["a", "b", "c", "d", "e", "f", "g"]
    ],
    "166": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "c", "e", "d", "g", "f", "h", "i"]
    ],
    "167": [
        ["a", "b", "c", "d", "e", "f"]
    ],
    "168": [
        ["a", "b", "c", "d"]
    ],
    "169": [
        ["a"]
    ],
    "170": [
        ["a"]
    ],
    "171": [
        ["a", "b", "c"]
    ],
    "172": [
        ["a", "b", "c"]
    ],
    "173": [
        ["a", "b", "c"]
    ],
    "174": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["a", "b", "e", "f", "c", "d", "g", "i", "h", "j", "k", "l"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "i", "k", "j", "l"],
        ["b", "a", "f", "e", "d", "c", "g", "i", "h", "k", "j", "l"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "i", "j", "k", "l"],
        ["c", "d", "e", "f", "a", "b", "h", "i", "g", "j", "k", "l"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i", "k", "j", "l"],
        ["d", "c", "f", "e", "b", "a", "h", "i", "g", "k", "j", "l"],
        ["e", "f", "a", "b", "c", "d", "i", "g", "h", "j", "k", "l"],
        ["e", "f", "c", "d", "a", "b", "i", "h", "g", "j", "k", "l"],
        ["f", "e", "b", "a", "d", "c", "i", "g", "h", "k", "j", "l"],
        ["f", "e", "d", "c", "b", "a", "i", "h", "g", "k", "j", "l"]
    ],
    "175": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["b", "a", "d", "c", "e", "g", "f", "h", "i", "k", "j", "l"]
    ],
    "176": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i"]
    ],
    "177": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"],
        ["b", "a", "d", "c", "e", "g", "f", "h", "i", "k", "j", "m", "l", "n"]
    ],
    "178": [
        ["a", "b", "c"]
    ],
    "179": [
        ["a", "b", "c"]
    ],
    "180": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["b", "a", "d", "c", "e", "f", "h", "g", "j", "i", "k"]
    ],
    "181": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["b", "a", "d", "c", "e", "f", "h", "g", "j", "i", "k"]
    ],
    "182": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i"]
    ],
    "183": [
        ["a", "b", "c", "d", "e", "f"]
    ],
    "184": [
        ["a", "b", "c", "d"]
    ],
    "185": [
        ["a", "b", "c", "d"]
    ],
    "186": [
        ["a", "b", "c", "d"]
    ],
    "187": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o"],
        ["a", "b", "e", "f", "c", "d", "g", "i", "h", "j", "k", "l", "m", "n", "o"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "i", "k", "j", "m", "l", "n", "o"],
        ["b", "a", "f", "e", "d", "c", "g", "i", "h", "k", "j", "m", "l", "n", "o"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "i", "j", "k", "l", "m", "n", "o"],
        ["c", "d", "e", "f", "a", "b", "h", "i", "g", "j", "k", "l", "m", "n", "o"],
        ["d", "c", "b", "a", "f", "e", "h", "g", "i", "k", "j", "m", "l", "n", "o"],
        ["d", "c", "f", "e", "b", "a", "h", "i", "g", "k", "j", "m", "l", "n", "o"],
        ["e", "f", "a", "b", "c", "d", "i", "g", "h", "j", "k", "l", "m", "n", "o"],
        ["e", "f", "c", "d", "a", "b", "i", "h", "g", "j", "k", "l", "m", "n", "o"],
        ["f", "e", "b", "a", "d", "c", "i", "g", "h", "k", "j", "m", "l", "n", "o"],
        ["f", "e", "d", "c", "b", "a", "i", "h", "g", "k", "j", "m", "l", "n", "o"]
    ],
    "188": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["a", "b", "e", "f", "c", "d", "g", "i", "h", "j", "k", "l"],
        ["c", "d", "a", "b", "e", "f", "h", "g", "i", "j", "k", "l"],
        ["c", "d", "e", "f", "a", "b", "h", "i", "g", "j", "k", "l"],
        ["e", "f", "a", "b", "c", "d", "i", "g", "h", "j", "k", "l"],
        ["e", "f", "c", "d", "a", "b", "i", "h", "g", "j", "k", "l"]
    ],
    "189": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["b", "a", "d", "c", "e", "g", "f", "h", "i", "k", "j", "l"]
    ],
    "190": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i"]
    ],
    "191": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r"],
        ["b", "a", "d", "c", "e", "g", "f", "h", "i", "k", "j", "m", "l", "n", "o", "q", "p", "r"]
    ],
    "192": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    ],
    "193": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]
    ],
    "194": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i", "j", "k", "l"]
    ],
    "195": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["a", "b", "c", "d", "e", "f", "h", "g", "i", "j"],
        ["b", "a", "d", "c", "e", "i", "g", "h", "f", "j"],
        ["b", "a", "d", "c", "e", "i", "h", "g", "f", "j"]
    ],
    "196": [
        ["a", "b", "c", "d", "e", "f", "g", "h"],
        ["a", "b", "d", "c", "e", "f", "g", "h"],
        ["b", "a", "c", "d", "e", "f", "g", "h"],
        ["b", "a", "d", "c", "e", "f", "g", "h"],
        ["c", "d", "a", "b", "e", "g", "f", "h"],
        ["c", "d", "b", "a", "e", "g", "f", "h"],
        ["d", "c", "a", "b", "e", "g", "f", "h"],
        ["d", "c", "b", "a", "e", "g", "f", "h"]
    ],
    "197": [
        ["a", "b", "c", "d", "e", "f"]
    ],
    "198": [
        ["a", "b"]
    ],
    "199": [
        ["a", "b", "c"]
    ],
    "200": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["a", "b", "c", "d", "e", "g", "f", "h", "i", "j", "k", "l"],
        ["b", "a", "d", "c", "h", "f", "g", "e", "i", "k", "j", "l"],
        ["b", "a", "d", "c", "h", "g", "f", "e", "i", "k", "j", "l"]
    ],
    "201": [
        ["a", "b", "c", "d", "e", "f", "g", "h"],
        ["a", "c", "b", "d", "e", "f", "g", "h"]
    ],
    "202": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "i"]
    ],
    "203": [
        ["a", "b", "c", "d", "e", "f", "g"],
        ["b", "a", "d", "c", "e", "f", "g"]
    ],
    "204": [
        ["a", "b", "c", "d", "e", "f", "g", "h"]
    ],
    "205": [
        ["a", "b", "c", "d"],
        ["b", "a", "c", "d"]
    ],
    "206": [
        ["a", "b", "c", "d", "e"],
        ["b", "a", "c", "d", "e"]
    ],
    "207": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "j", "i", "k"]
    ],
    "208": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"],
        ["a", "b", "c", "d", "f", "e", "g", "h", "j", "i", "k", "l", "m"],
        ["a", "c", "b", "d", "e", "f", "g", "h", "i", "j", "l", "k", "m"],
        ["a", "c", "b", "d", "f", "e", "g", "h", "j", "i", "l", "k", "m"]
    ],
    "209": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["b", "a", "c", "d", "e", "f", "h", "g", "i", "j"]
    ],
    "210": [
        ["a", "b", "c", "d", "e", "f", "g", "h"],
        ["b", "a", "d", "c", "e", "f", "g", "h"]
    ],
    "211": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]
    ],
    "212": [
        ["a", "b", "c", "d", "e"],
        ["b", "a", "c", "d", "e"]
    ],
    "213": [
        ["a", "b", "c", "d", "e"],
        ["b", "a", "c", "d", "e"]
    ],
    "214": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "d", "c", "e", "f", "h", "g", "i"]
    ],
    "215": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"],
        ["b", "a", "d", "c", "e", "g", "f", "h", "i", "j"]
    ],
    "216": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["a", "b", "d", "c", "e", "f", "g", "h", "i"],
        ["b", "a", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "i"],
        ["c", "d", "a", "b", "e", "g", "f", "h", "i"],
        ["c", "d", "b", "a", "e", "g", "f", "h", "i"],
        ["d", "c", "a", "b", "e", "g", "f", "h", "i"],
        ["d", "c", "b", "a", "e", "g", "f", "h", "i"]
    ],
    "217": [
        ["a", "b", "c", "d", "e", "f", "g", "h"]
    ],
    "218": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["a", "b", "d", "c", "e", "f", "h", "g", "i"]
    ],
    "219": [
        ["a", "b", "c", "d", "e", "f", "g", "h"],
        ["b", "a", "d", "c", "e", "g", "f", "h"]
    ],
    "220": [
        ["a", "b", "c", "d", "e"],
        ["b", "a", "c", "d", "e"]
    ],
    "221": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"],
        ["b", "a", "d", "c", "f", "e", "g", "h", "j", "i", "l", "k", "m", "n"]
    ],
    "222": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
    ],
    "223": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["a", "b", "d", "c", "e", "f", "h", "g", "i", "j", "k", "l"]
    ],
    "224": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["a", "c", "b", "d", "e", "f", "g", "h", "j", "i", "k", "l"]
    ],
    "225": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
        ["b", "a", "c", "d", "e", "f", "g", "i", "h", "j", "k", "l"]
    ],
    "226": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]
    ],
    "227": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        ["b", "a", "d", "c", "e", "f", "g", "h", "i"]
    ],
    "228": [
        ["a", "b", "c", "d", "e", "f", "g", "h"]
    ],
    "229": [
        ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]
    ],
    "230": [
        ["a", "b", "c", "d", "e", "f", "g", "h"]
    ]
}

SPACE_GROUPS_FOR_EACH_PEARSON = {
    "aP":[1,2],
    "mP":[3,4,6,7,10,11,13,14],
    "mC":[5,8,9,12,15],
    "oP":[6,17,18,19,25,26,27,28,29,30,31,32,33,34,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62],
    "oC":[20,21,35,36,37,38,39,40,41,63,64,65,66,67,68],
    "oI":[23,24,44,45,46,71,72,73,74],
    "oF":[22,42,43,69,70],
    "tP":[75,76,77,78,81,83,84,85,86,89,90,91,92,93,94,95,96,99,100,101,102,103,104,105,106,111,112,113,114,115,116,117,118,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138],
    "tI":[79,80,82,87,88,97,98,107,108,109,110,119,120,121,122,139,140,141,142],
    "hP":[143,144,145,147,149,150,151,152,153,154,156,157,158,159,162,163,164,165,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194],
    "hR":[146,148,155,160,161,166,167],
    "cP":[195,198,200,201,205,207,208,212,213,215,218,221,222,223,224],
    "cF":[196,202,203,209,210,216,219,225,226,227,228],
    "cI":[197,199,204,206,211,214,217,220,229,230]
}

def get_pearson_from_space_group(sgnum:Union[int,str]):
    for pearson in SPACE_GROUPS_FOR_EACH_PEARSON:
        if int(sgnum) in SPACE_GROUPS_FOR_EACH_PEARSON[pearson]:
            return pearson
        
def get_formal_pearson_from_space_group(sgnum:Union[int,str]):
    """
    same as above, except distinguish between "oA" and "oC"
    """
    pearson = get_pearson_from_space_group(sgnum)
    if pearson == "oC":
        if sgnum in A_CENTERED_ORTHORHOMBIC_GROUPS:
            return "oA"
    return pearson

POSSIBLE_PRIMITIVE_SHIFTS = {
    1: [
        [0,0,0]
    ],
    2: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    3: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0,0],
        [0.5,0,0.5]
    ],
    4: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    5: [
        [0,0,0],
        [0,0,0.5]
    ],
    6: [
        [0,0,0],
        [0,0.5,0]
    ],
    7: [
        [0,0,0],
        [0,0.5,0],
        [0,0,0.5],
        [0,0.5,0.5]
    ],
    8: [
        [0,0,0]
    ],
    9: [
        [0,0,0],
        [0.75,0.25,0],
        [0.5,0.5,0],
        [0.25,0.75,0],
        [0,0,0.5],
        [0.75,0.25,0.5],
        [0.5,0.5,0.5],
        [0.25,0.75,0.5]
    ],
    10: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    11: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    12: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    13: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    14: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    15: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5]
    ],
    16: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    17: [
        [0,0,0],
        [0,0,0.25],
        [0,0,0.5],
        [0,0,0.75],
        [0,0.5,0],
        [0,0.5,0.25],
        [0,0.5,0.5],
        [0,0.5,0.75],
        [0.5,0,0],
        [0.5,0,0.25],
        [0.5,0,0.5],
        [0.5,0,0.75],
        [0.5,0.5,0],
        [0.5,0.5,0.25],
        [0.5,0.5,0.5],
        [0.5,0.5,0.75]
    ],
    18: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    19: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.25,0.25,0.25],
        [0.25,0.25,0.75],
        [0.25,0.75,0.25],
        [0.25,0.75,0.75],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0.75,0.25,0.25],
        [0.75,0.25,0.75],
        [0.75,0.75,0.25],
        [0.75,0.75,0.75]
    ],
    20: [
        [0,0,0],
        [0,0,0.25],
        [0,0,0.5],
        [0,0,0.75],
        [0.5,0.5,0],
        [0.5,0.5,0.25],
        [0.5,0.5,0.5],
        [0.5,0.5,0.75]
    ],
    21: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    22: [
        [0,0,0],
        [0.25,0.25,0.25],
        [0.5,0.5,0.5],
        [0.75,0.75,0.75]
    ],
    23: [
        [0,0,0],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [0.5,0.5,0]
    ],
    24: [
        [0,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0.5,0],
        [0.5,0,0],
        [0,0,0.5]
    ],
    25: [
        [0,0,0],
        [0,0.5,0],
        [0.5,0,0],
        [0.5,0.5,0]
    ],
    26: [
        [0,0,0],
        [0,0.5,0],
        [0.5,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0,0.5,0.5],
        [0.5,0,0.5],
        [0.5,0.5,0.5]
    ],
    27: [
        [0,0,0],
        [0,0.5,0],
        [0.5,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0,0.5,0.5],
        [0.5,0,0.5],
        [0.5,0.5,0.5]
    ],
    28: [
        [0,0,0],
        [0,0.5,0],
        [0.5,0,0],
        [0.5,0.5,0]
    ],
    29: [
        [0,0,0],
        [0,0.5,0],
        [0.5,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0,0.5,0.5],
        [0.5,0,0.5],
        [0.5,0.5,0.5]
    ],
    30: [
        [0,0,0],
        [0,0.5,0],
        [0.5,0,0],
        [0.5,0.5,0],
        [0,0.5,0.5],
        [0,0,0.5],
        [0.5,0.5,0.5],
        [0.5,0,0.5]
    ],
    31: [
        [0,0,0],
        [0,0.5,0],
        [0.5,0,0],
        [0.5,0.5,0],
        [0.5,0,0.5],
        [0.5,0.5,0.5],
        [0,0,0.5],
        [0,0.5,0.5]
    ],
    32: [
        [0,0,0],
        [0,0.5,0],
        [0.5,0,0],
        [0.5,0.5,0]
    ],
    33: [
        [0,0,0],
        [0,0.5,0],
        [0.5,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0,0.5,0.5],
        [0.5,0,0.5],
        [0.5,0.5,0.5]
    ],
    34: [
        [0,0,0],
        [0,0.5,0],
        [0.5,0,0],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [0,0,0.5]
    ],
    35: [
        [0,0,0],
        [0.5,0.5,0]
    ],
    36: [
        [0,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.5,0.5,0.5]
    ],
    37: [
        [0,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.5,0.5,0.5]
    ],
    38: [
        [0,0,0],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0,0.5]
    ],
    39: [
        [0,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.5,0.5,0.5]
    ],
    40: [
        [0,0,0],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0,0.5]
    ],
    41: [
        [0,0,0],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0,0.5]
    ],
    42: [
        [0,0,0]
    ],
    43: [
        [0,0,0],
        [0,0,0.5],
        [0.25,0.25,0.25],
        [0.25,0.25,0.75]
    ],
    44: [
        [0,0,0],
        [0,0.5,0.5]
    ],
    45: [
        [0,0,0],
        [0,0.5,0.5],
        [0.5,0.5,0],
        [0.5,0,0.5]
    ],
    46: [
        [0,0,0],
        [0,0.5,0.5]
    ],
    47: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    48: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    49: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    50: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    51: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    52: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    53: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    54: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    55: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    56: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    57: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    58: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    59: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    60: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    61: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    62: [
        [0,0,0],
        [0,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    63: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    64: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    65: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    66: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    67: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0,0.5,0.5]
    ],
    68: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0,0.5,0.5]
    ],
    69: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    70: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0.5,0,0],
        [0,0.5,0.5],
        [0,0.5,0],
        [0.5,0,0.5],
        [0,0,0.5],
        [0.5,0.5,0]
    ],
    71: [
        [0,0,0],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [0.5,0.5,0]
    ],
    72: [
        [0,0,0],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [0.5,0.5,0]
    ],
    73: [
        [0,0,0],
        [0.5,0,0.5],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0.5,0],
        [0.5,0,0],
        [0,0,0.5]
    ],
    74: [
        [0,0,0],
        [0.5,0,0.5],
        [0,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0.5],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.5,0,0]
    ],
    75: [
        [0,0,0],
        [0.5,0.5,0]
    ],
    76: [
        [0,0,0],
        [0.5,0.5,0],
        [0,0,0.25],
        [0.5,0.5,0.25],
        [0,0,0.5],
        [0.5,0.5,0.5],
        [0,0,0.75],
        [0.5,0.5,0.75]
    ],
    77: [
        [0,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.5,0.5,0.5]
    ],
    78: [
        [0,0,0],
        [0.5,0.5,0],
        [0,0,0.25],
        [0.5,0.5,0.25],
        [0,0,0.5],
        [0.5,0.5,0.5],
        [0,0,0.75],
        [0.5,0.5,0.75]
    ],
    79: [
        [0,0,0]
    ],
    80: [
        [0,0,0],
        [0,0.5,0.5],
        [0.75,0.25,0.5],
        [0.75,0.75,0]
    ],
    81: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    82: [
        [0,0,0],
        [0.5,0.5,0],
        [0.25,0.75,0.5],
        [0.75,0.25,0.5]
    ],
    83: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    84: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    85: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5]
    ],
    86: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0.5],
        [0,0.5,0],
        [0.5,0,0.5],
        [0.5,0,0]
    ],
    87: [
        [0,0,0],
        [0.5,0.5,0]
    ],
    88: [
        [0,0,0],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0,0.5],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [0,0.5,0],
        [0.5,0,0]
    ],
    89: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    90: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    91: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0,0.25],
        [0,0,0.75],
        [0.5,0.5,0.25],
        [0.5,0.5,0.75]
    ],
    92: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0.5,0.5,0.25],
        [0.5,0.5,0.75],
        [0,0,0.25],
        [0,0,0.75]
    ],
    93: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    94: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    95: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0,0.25],
        [0,0,0.75],
        [0.5,0.5,0.25],
        [0.5,0.5,0.75]
    ],
    96: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0.5,0.5,0.25],
        [0.5,0.5,0.75],
        [0,0,0.25],
        [0,0,0.75]
    ],
    97: [
        [0,0,0],
        [0.5,0.5,0]
    ],
    98: [
        [0,0,0],
        [0.5,0.5,0],
        [0.25,0.75,0.5],
        [0.75,0.25,0.5]
    ],
    99: [
        [0,0,0],
        [0.5,0.5,0]
    ],
    100: [
        [0,0,0],
        [0.5,0.5,0]
    ],
    101: [
        [0,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.5,0.5,0.5]
    ],
    102: [
        [0,0,0],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0,0.5]
    ],
    103: [
        [0,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.5,0.5,0.5]
    ],
    104: [
        [0,0,0],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0,0.5]
    ],
    105: [
        [0,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.5,0.5,0.5]
    ],
    106: [
        [0,0,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.5,0.5,0.5]
    ],
    107: [
        [0,0,0]
    ],
    108: [
        [0,0,0],
        [0.5,0.5,0]
    ],
    109: [
        [0,0,0],
        [0,0.5,0.5],
        [0.75,0.25,0.5],
        [0.75,0.75,0]
    ],
    110: [
        [0,0,0],
        [0,0.5,0.5],
        [0.5,0.5,0],
        [0.5,0,0.5],
        [0.75,0.25,0.5],
        [0.75,0.75,0],
        [0.25,0.75,0.5],
        [0.25,0.25,0]
    ],
    111: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    112: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    113: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    114: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    115: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    116: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    117: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    118: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    119: [
        [0,0,0],
        [0.5,0.5,0],
        [0.25,0.75,0.5],
        [0.75,0.25,0.5]
    ],
    120: [
        [0,0,0],
        [0.5,0.5,0],
        [0.25,0.75,0.5],
        [0.75,0.25,0.5]
    ],
    121: [
        [0,0,0],
        [0.5,0.5,0]
    ],
    122: [
        [0,0,0],
        [0.5,0.5,0],
        [0.25,0.75,0.5],
        [0.75,0.25,0.5]
    ],
    123: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    124: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    125: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5]
    ],
    126: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5]
    ],
    127: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    128: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    129: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5]
    ],
    130: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5]
    ],
    131: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    132: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    133: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5]
    ],
    134: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0.5],
        [0,0.5,0],
        [0.5,0,0.5],
        [0.5,0,0]
    ],
    135: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    136: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5]
    ],
    137: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5]
    ],
    138: [
        [0,0,0],
        [0,0,0.5],
        [0.5,0.5,0],
        [0.5,0.5,0.5],
        [0,0.5,0.5],
        [0,0.5,0],
        [0.5,0,0.5],
        [0.5,0,0]
    ],
    139: [
        [0,0,0],
        [0.5,0.5,0]
    ],
    140: [
        [0,0,0],
        [0.5,0.5,0]
    ],
    141: [
        [0,0,0],
        [0.5,0.5,0],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [0,0,0.5],
        [0.5,0.5,0.5],
        [0,0.5,0],
        [0.5,0,0]
    ],
    142: [
        [0,0,0],
        [0.5,0.5,0],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [0.5,0.5,0.5],
        [0,0,0.5],
        [0,0.5,0],
        [0.5,0,0]
    ],
    143: [
        [0,0,0],
        [0.3333333333333333,0.6666666666666666,0],
        [0.6666666666666666,0.3333333333333333,0]
    ],
    144: [
        [0,0,0],
        [0.3333333333333333,0.6666666666666666,0],
        [0.6666666666666666,0.3333333333333333,0],
        [0,0,0.3333333333333333],
        [0.3333333333333333,0.6666666666666666,0.3333333333333333],
        [0.6666666666666666,0.3333333333333333,0.3333333333333333],
        [0,0,0.6666666666666666],
        [0.3333333333333333,0.6666666666666666,0.6666666666666666],
        [0.6666666666666666,0.3333333333333333,0.6666666666666666]
    ],
    145: [
        [0,0,0],
        [0.3333333333333333,0.6666666666666666,0],
        [0.6666666666666666,0.3333333333333333,0],
        [0,0,0.3333333333333333],
        [0.3333333333333333,0.6666666666666666,0.3333333333333333],
        [0.6666666666666666,0.3333333333333333,0.3333333333333333],
        [0,0,0.6666666666666666],
        [0.3333333333333333,0.6666666666666666,0.6666666666666666],
        [0.6666666666666666,0.3333333333333333,0.6666666666666666]
    ],
    146: [
        [0,0,0]
    ],
    147: [
        [0,0,0],
        [0,0,0.5]
    ],
    148: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    149: [
        [0,0,0],
        [0,0,0.5],
        [0.3333333333333333,0.6666666666666666,0],
        [0.3333333333333333,0.6666666666666666,0.5],
        [0.6666666666666666,0.3333333333333333,0],
        [0.6666666666666666,0.3333333333333333,0.5]
    ],
    150: [
        [0,0,0],
        [0,0,0.5]
    ],
    151: [
        [0,0,0],
        [0,0,0.5],
        [0.3333333333333333,0.6666666666666666,0],
        [0.3333333333333333,0.6666666666666666,0.5],
        [0.6666666666666666,0.3333333333333333,0],
        [0.6666666666666666,0.3333333333333333,0.5],
        [0,0,0.3333333333333333],
        [0,0,0.8333333333333333],
        [0.3333333333333333,0.6666666666666666,0.3333333333333333],
        [0.3333333333333333,0.6666666666666666,0.8333333333333333],
        [0.6666666666666666,0.3333333333333333,0.3333333333333333],
        [0.6666666666666666,0.3333333333333333,0.8333333333333333],
        [0,0,0.6666666666666666],
        [0,0,0.16666666666666652],
        [0.3333333333333333,0.6666666666666666,0.6666666666666666],
        [0.3333333333333333,0.6666666666666666,0.16666666666666652],
        [0.6666666666666666,0.3333333333333333,0.6666666666666666],
        [0.6666666666666666,0.3333333333333333,0.16666666666666652]
    ],
    152: [
        [0,0,0],
        [0,0,0.5],
        [0,0,0.3333333333333333],
        [0,0,0.8333333333333333],
        [0,0,0.6666666666666666],
        [0,0,0.16666666666666652]
    ],
    153: [
        [0,0,0],
        [0,0,0.5],
        [0.3333333333333333,0.6666666666666666,0],
        [0.3333333333333333,0.6666666666666666,0.5],
        [0.6666666666666666,0.3333333333333333,0],
        [0.6666666666666666,0.3333333333333333,0.5],
        [0,0,0.3333333333333333],
        [0,0,0.8333333333333333],
        [0.3333333333333333,0.6666666666666666,0.3333333333333333],
        [0.3333333333333333,0.6666666666666666,0.8333333333333333],
        [0.6666666666666666,0.3333333333333333,0.3333333333333333],
        [0.6666666666666666,0.3333333333333333,0.8333333333333333],
        [0,0,0.6666666666666666],
        [0,0,0.16666666666666652],
        [0.3333333333333333,0.6666666666666666,0.6666666666666666],
        [0.3333333333333333,0.6666666666666666,0.16666666666666652],
        [0.6666666666666666,0.3333333333333333,0.6666666666666666],
        [0.6666666666666666,0.3333333333333333,0.16666666666666652]
    ],
    154: [
        [0,0,0],
        [0,0,0.5],
        [0,0,0.3333333333333333],
        [0,0,0.8333333333333333],
        [0,0,0.6666666666666666],
        [0,0,0.16666666666666652]
    ],
    155: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    156: [
        [0,0,0],
        [0.3333333333333333,0.6666666666666666,0],
        [0.6666666666666666,0.3333333333333333,0]
    ],
    157: [
        [0,0,0]
    ],
    158: [
        [0,0,0],
        [0.3333333333333333,0.6666666666666666,0],
        [0.6666666666666666,0.3333333333333333,0],
        [0,0,0.5],
        [0.3333333333333333,0.6666666666666666,0.5],
        [0.6666666666666666,0.3333333333333333,0.5]
    ],
    159: [
        [0,0,0],
        [0,0,0.5]
    ],
    160: [
        [0,0,0]
    ],
    161: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    162: [
        [0,0,0],
        [0,0,0.5]
    ],
    163: [
        [0,0,0],
        [0,0,0.5]
    ],
    164: [
        [0,0,0],
        [0,0,0.5]
    ],
    165: [
        [0,0,0],
        [0,0,0.5]
    ],
    166: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    167: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    168: [
        [0,0,0]
    ],
    169: [
        [0,0,0],
        [0,0,0.16666666666666666],
        [0,0,0.3333333333333333],
        [0,0,0.5],
        [0,0,0.6666666666666666],
        [0,0,0.8333333333333334]
    ],
    170: [
        [0,0,0],
        [0,0,0.16666666666666666],
        [0,0,0.3333333333333333],
        [0,0,0.5],
        [0,0,0.6666666666666666],
        [0,0,0.8333333333333334]
    ],
    171: [
        [0,0,0],
        [0,0,0.3333333333333333],
        [0,0,0.6666666666666666]
    ],
    172: [
        [0,0,0],
        [0,0,0.3333333333333333],
        [0,0,0.6666666666666666]
    ],
    173: [
        [0,0,0],
        [0,0,0.5]
    ],
    174: [
        [0,0,0],
        [0,0,0.5],
        [0.3333333333333333,0.6666666666666666,0],
        [0.3333333333333333,0.6666666666666666,0.5],
        [0.6666666666666666,0.3333333333333333,0],
        [0.6666666666666666,0.3333333333333333,0.5]
    ],
    175: [
        [0,0,0],
        [0,0,0.5]
    ],
    176: [
        [0,0,0],
        [0,0,0.5]
    ],
    177: [
        [0,0,0],
        [0,0,0.5]
    ],
    178: [
        [0,0,0],
        [0,0,0.5],
        [0,0,0.16666666666666666],
        [0,0,0.6666666666666666],
        [0,0,0.3333333333333333],
        [0,0,0.8333333333333333]
    ],
    179: [
        [0,0,0],
        [0,0,0.5],
        [0,0,0.16666666666666666],
        [0,0,0.6666666666666666],
        [0,0,0.3333333333333333],
        [0,0,0.8333333333333333]
    ],
    180: [
        [0,0,0],
        [0,0,0.5],
        [0,0,0.3333333333333333],
        [0,0,0.8333333333333333],
        [0,0,0.6666666666666666],
        [0,0,0.16666666666666652]
    ],
    181: [
        [0,0,0],
        [0,0,0.5],
        [0,0,0.3333333333333333],
        [0,0,0.8333333333333333],
        [0,0,0.6666666666666666],
        [0,0,0.16666666666666652]
    ],
    182: [
        [0,0,0],
        [0,0,0.5]
    ],
    183: [
        [0,0,0]
    ],
    184: [
        [0,0,0],
        [0,0,0.5]
    ],
    185: [
        [0,0,0],
        [0,0,0.5]
    ],
    186: [
        [0,0,0],
        [0,0,0.5]
    ],
    187: [
        [0,0,0],
        [0,0,0.5],
        [0.3333333333333333,0.6666666666666666,0],
        [0.3333333333333333,0.6666666666666666,0.5],
        [0.6666666666666666,0.3333333333333333,0],
        [0.6666666666666666,0.3333333333333333,0.5]
    ],
    188: [
        [0,0,0],
        [0,0,0.5],
        [0.3333333333333333,0.6666666666666666,0],
        [0.3333333333333333,0.6666666666666666,0.5],
        [0.6666666666666666,0.3333333333333333,0],
        [0.6666666666666666,0.3333333333333333,0.5]
    ],
    189: [
        [0,0,0],
        [0,0,0.5]
    ],
    190: [
        [0,0,0],
        [0,0,0.5]
    ],
    191: [
        [0,0,0],
        [0,0,0.5]
    ],
    192: [
        [0,0,0],
        [0,0,0.5]
    ],
    193: [
        [0,0,0],
        [0,0,0.5]
    ],
    194: [
        [0,0,0],
        [0,0,0.5]
    ],
    195: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    196: [
        [0,0,0],
        [0.25,0.25,0.25],
        [0.5,0.5,0.5],
        [0.75,0.75,0.75]
    ],
    197: [
        [0,0,0]
    ],
    198: [
        [0,0,0],
        [0.25,0.25,0.25],
        [0.5,0.5,0.5],
        [0.75,0.75,0.75],
        [0,0.5,0.5],
        [0.25,0.75,0.75],
        [0.5,0,0],
        [0.75,0.25,0.25],
        [0.5,0,0.5],
        [0.75,0.25,0.75],
        [0,0.5,0],
        [0.25,0.75,0.25],
        [0.5,0.5,0],
        [0.75,0.75,0.25],
        [0,0,0.5],
        [0.25,0.25,0.75]
    ],
    199: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.5,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0]
    ],
    200: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    201: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0,0.5,0],
        [0.5,0.5,0],
        [0,0,0.5]
    ],
    202: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    203: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0.5,0,0],
        [0,0.5,0.5],
        [0,0.5,0],
        [0.5,0,0.5],
        [0,0,0.5],
        [0.5,0.5,0]
    ],
    204: [
        [0,0,0]
    ],
    205: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0,0.5,0],
        [0.5,0.5,0],
        [0,0,0.5]
    ],
    206: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.5,0,0.5],
        [0,0.5,0],
        [0,0.5,0.5],
        [0.5,0,0]
    ],
    207: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    208: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    209: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    210: [
        [0,0,0],
        [0.25,0.25,0.25],
        [0.5,0.5,0.5],
        [0.75,0.75,0.75]
    ],
    211: [
        [0,0,0]
    ],
    212: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.25,0.25,0.25],
        [0.75,0.75,0.75],
        [0.25,0.75,0.75],
        [0.75,0.25,0.25],
        [0.5,0,0.5],
        [0,0.5,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.75,0.25,0.75],
        [0.25,0.75,0.25],
        [0.75,0.75,0.25],
        [0.25,0.25,0.75]
    ],
    213: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.25,0.25,0.75],
        [0.75,0.75,0.25],
        [0.25,0.75,0.25],
        [0.75,0.25,0.75],
        [0.5,0,0.5],
        [0,0.5,0],
        [0.5,0.5,0],
        [0,0,0.5],
        [0.75,0.25,0.25],
        [0.25,0.75,0.75],
        [0.75,0.75,0.75],
        [0.25,0.25,0.25]
    ],
    214: [
        [0,0,0],
        [0.5,0.5,0],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [0.5,0.5,0.5],
        [0,0,0.5],
        [0,0.5,0],
        [0.5,0,0]
    ],
    215: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    216: [
        [0,0,0],
        [0.25,0.25,0.25],
        [0.5,0.5,0.5],
        [0.75,0.75,0.75]
    ],
    217: [
        [0,0,0]
    ],
    218: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    219: [
        [0,0,0],
        [0.25,0.25,0.25],
        [0.5,0.5,0.5],
        [0.75,0.75,0.75]
    ],
    220: [
        [0,0,0],
        [0.5,0.5,0],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [0.5,0.5,0.5],
        [0,0,0.5],
        [0,0.5,0],
        [0.5,0,0]
    ],
    221: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    222: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0,0,0.5],
        [0.5,0.5,0],
        [0,0.5,0],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [0.5,0,0]
    ],
    223: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    224: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0,0.5,0.5],
        [0.5,0,0],
        [0.5,0,0.5],
        [0,0.5,0],
        [0.5,0.5,0],
        [0,0,0.5]
    ],
    225: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    226: [
        [0,0,0],
        [0.5,0.5,0.5]
    ],
    227: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0.5,0,0],
        [0,0.5,0.5],
        [0,0.5,0],
        [0.5,0,0.5],
        [0,0,0.5],
        [0.5,0.5,0]
    ],
    228: [
        [0,0,0],
        [0.5,0.5,0.5],
        [0.5,0,0],
        [0,0.5,0.5],
        [0,0.5,0],
        [0.5,0,0.5],
        [0,0,0.5],
        [0.5,0.5,0]
    ],
    229: [
        [0,0,0]
    ],
    230: [
        [0,0,0],
        [0.5,0.5,0],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [0.5,0.5,0.5],
        [0,0,0.5],
        [0,0.5,0],
        [0.5,0,0]
    ]
}
