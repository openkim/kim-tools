"""Crystal Symmetry utilities that implement functionality not present in AFLOW"""

from numpy.typing import ArrayLike
from typing import Optional,List
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
]

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
    