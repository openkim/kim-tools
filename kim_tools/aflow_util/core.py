"""Tools for working with crystal prototypes using the AFLOW command line tool"""
import numpy as np
from numpy.typing import ArrayLike
import json
import subprocess
import sys
import os
import ase
from ase.cell import Cell
from ase import Atoms
import ase.spacegroup
from ase.spacegroup.symmetrize import check_symmetry
from curses.ascii import isalpha, isupper, isdigit
from typing import Dict, List, Tuple, Union, Optional
from tempfile import NamedTemporaryFile
from sympy import parse_expr,matrix2numpy,linear_eq_to_matrix,Symbol
from dataclasses import dataclass
from ..symmetry_util import are_in_same_wyckoff_set, WYCKOFF_MULTIPLICITIES, CENTERING_DIVISORS
from operator import attrgetter
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename='kim-tools.log',level=logging.INFO,force=True)


__author__ = ["ilia Nikiforov", "Ellad Tadmor"]
__all__ = [
    "incorrectNumAtomsException",
    "EquivalentEqnSet",
    "EquivalentAtomSet",
    "group_positions_by_wyckoff",
    "spglibFailureException",
    "incorrectSpaceGroupException",
    "incorrectNumSpeciesException",
    "inconsistentWyckoffException",    
    "CENTERING_DIVISORS",
    "check_number_of_atoms",
    "split_parameter_array",
    "internal_parameter_sort_key",
    "get_stoich_reduced_list_from_prototype",
    "get_wyckoff_lists_from_prototype",
    "get_space_group_number_from_prototype",
    "get_pearson_symbol_from_prototype",
    "get_centering_divisor_from_prototype",
    "read_shortnames",
    "get_real_to_virtual_species_map",
    "solve_for_free_params",
    "AFLOW"
]

class incorrectNumAtomsException(Exception):
    """
    Raised when the number of atoms in the atoms object or in the WYCCAR returned by aflow --sgdata does not match the number of atoms in the Pearson symbol of the prototype label
    """

class spglibFailureException(Exception):
    """
    Raised when an issue with the spglib analysis occurs
    """

class incorrectSpaceGroupException(Exception):
    """
    Raised when spglib or aflow --sgdata detects a different space group than the one specified in the prototype label
    """
    
class incorrectNumSpeciesException(Exception):
    """
    Raised when number of species is inconsistent
    """
    
class inconsistentWyckoffException(Exception):
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
    param_names: List[str] # The n free parameters associated with this Wyckoff postition, 0 <= n <= 3
    coeff_matrix_list: List[ArrayLike] # m x 3 x n matrices of coefficients, where m is the multiplicity of the Wyckoff position
    const_terms_list: List[ArrayLike] # m x 3 x 1 columns of constant terms in the coordinates. This gets subtracted from the RHS when solving
    

@dataclass
class EquivalentAtomSet:
    """
    Set of equivalent atoms
    """
    species: str
    wyckoff_letter: str
    frac_position_list: List[ArrayLike] # m x 3 x 1 columns
    
def check_number_of_atoms(atoms: Atoms, prototype_label: str, primitive_cell: bool = True) -> None:
    """
    Check if the Atoms object (which must be a conventional or primitive unit cell)
    has the correct number of atoms according to prototype_label
    
    Raises: 
    incorrectNumAtomsException
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
    
    if centering == 'R':
        num_conv_cell *= 3

    if (not primitive_cell):
        num_lattice = 1
    else:
        num_lattice = CENTERING_DIVISORS[centering]
        
    # This check is probably really extraneous, but better safe than sorry
    if num_conv_cell % num_lattice != 0:
        raise incorrectNumAtomsException("WARNING: Number of atoms in conventional cell %d derived from Pearson symbol of prototype %s is not divisible by the number of lattice points %d"%(num_conv_cell,prototype_label,num_lattice))
        
    num_cell = num_conv_cell/num_lattice
            
    if len(atoms)!=num_cell:
        raise incorrectNumAtomsException("WARNING: Number of ASE atoms %d does not match Pearson symbol of prototype %s"%(len(atoms),prototype_label))

def split_parameter_array(parameter_names: List[str], list_to_split: Optional[List] = None) -> Tuple[List,List]:
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
    
    assert (len(list_to_split) == len(parameter_names)), \
        "`list_to_split` must have the same length as `parameter_names`"
    
    in_internal_part = False
    
    cell_part = []
    internal_part = []
    
    CARTESIAN_AXES = ['x','y','z']
    
    for name,value in zip(parameter_names,list_to_split):
        assert (isinstance(name,str)), "At least one element of `parameter_names` is not a string."
        if not in_internal_part: 
            if name[0] in CARTESIAN_AXES:
                in_internal_part = True
        else: # means we have already encountered an internal coordinate in a past iteration
            assert (name[0] in CARTESIAN_AXES), \
                "`parameter_names` seems to have an internal parameter followed by a non-internal one"
        
        if in_internal_part:
            internal_part.append(value)
        else:
            cell_part.append(value)
            
    return cell_part, internal_part

def internal_parameter_sort_key(parameter_name: Union[Symbol,str] ) -> int:
    """
    Sorting key for internal free parameters. Sort by number first, then letter
    """
    parameter_name_str = str(parameter_name)
    axis = parameter_name_str[0]
    assert axis == 'x' or axis == 'y' or axis == 'z', 'Parameter name must start with x, y, or z'
    number = int(parameter_name_str[1:])
    return 1000*number + ord(axis)

def group_positions_by_wyckoff(atoms: Atoms, prototype_label: Optional[str] = None) -> List[EquivalentAtomSet]:
    """
    Return a list of objects representing sets of equivalent atoms
    
    Args:
        atoms: atoms object to analyze
        prototype_label:
            If this is provided, a consistency check will be made.
            It is assumed that the virtual species in the prototype label
            are alphabetized commensurately with the real species in ``atoms``,
            and that the Wyckoff positions within each species are alphabetized
            as well (as it should be in a valid prototype label)
            TODO: Make a checker for "valid prototype label"
            
    Raises:
        spglibFailureException: if spglib fails
        incorrectSpaceGroupException: if ``prototype_label`` is provided and disagrees with spglib space group
        incorrectNumSpeciesException: if ``prototype_label`` is provided and disagrees with number of species in ``atoms``
        inconsistentWyckoffException: if ``prototype_label`` is provided and disagrees with spglib Wyckoffs
    Returns:
        List of object each containing string attributes ``species``, ``wyckoff_letter``,
        and list of 3x1 arrays ``frac_position_list``
    """
    if prototype_label is not None:
        check_number_of_atoms(atoms,prototype_label)
    
    try:
        ds = check_symmetry(atoms)
    except Exception as e:
        raise spglibFailureException(f'spglib encountered the following exception:\n{e}')
    
    if ds is None:
        raise spglibFailureException(f'spglib returned ``None`` dataset')
    
    inequivalent_atom_indices = sorted(list(set(ds.crystallographic_orbits)))
    
    # initialize return list    
    equivalent_atom_set_list = []    
    for inequivalent_atom_index in inequivalent_atom_indices:
        equivalent_atom_set_list.append(
            EquivalentAtomSet(
                atoms.get_chemical_symbols()[inequivalent_atom_index],
                ds.wyckoffs[inequivalent_atom_index],
                []
            )
        )
    
    # fill with coordinates
    for frac_pos,repr_atom_index in zip(atoms.get_scaled_positions(),ds.crystallographic_orbits):
        for i,equivalent_atom_set in enumerate(equivalent_atom_set_list):
            if inequivalent_atom_indices[i] == repr_atom_index:
                equivalent_atom_set.frac_position_list.append(frac_pos.reshape(3,1))
                continue

    # sort by letter first then species
    equivalent_atom_set_list.sort(key=attrgetter('wyckoff_letter'))
    equivalent_atom_set_list.sort(key=attrgetter('species'))
    
    # check consistency with prototype label
    if prototype_label is not None:
        space_group_number = get_space_group_number_from_prototype(prototype_label)
        if ds.number != space_group_number:
            raise incorrectSpaceGroupException(
                f'spglib detected space group {ds.number}, against label {space_group_number}')
        wyckoff_lists = get_wyckoff_lists_from_prototype(prototype_label)
        if len(wyckoff_lists) != len(set(atoms.get_chemical_symbols())):
            raise incorrectNumSpeciesException(
                f'Prototype label {prototype_label} has {len(wyckoff_lists)} species-Wyckoff sections\n'
                f'but ``atoms`` has {len(set(atoms.get_chemical_symbols()))} unique species')
        wyckoff_lists_concatenated = ''
        for wyckoff_list in wyckoff_lists:
            wyckoff_lists_concatenated += wyckoff_list    
        if len(wyckoff_lists_concatenated) != len(equivalent_atom_set_list):
            raise inconsistentWyckoffException(
                f'Prototype label {prototype_label} indicates {len(wyckoff_lists_concatenated)} '
                f'inequivalent atoms but I found {len(equivalent_atom_set_list)}')
        # everything should be ordered consistently with each other,
        # except within a Wyckoff set which we will have to map later
        for label_wyckoff_letter,equivalent_atom_set in \
            zip(wyckoff_lists_concatenated,equivalent_atom_set_list):
                if not are_in_same_wyckoff_set(equivalent_atom_set.wyckoff_letter,label_wyckoff_letter,space_group_number):
                    raise inconsistentWyckoffException(
                        f'Prototype Wyckoff letter {label_wyckoff_letter} and '
                        f'spglib Wyckoff letter {equivalent_atom_set.wyckoff_letter} '
                        'are not in the same Wyckoff set'
                    )
                if label_wyckoff_letter != equivalent_atom_set.wyckoff_letter:
                    logger.info(f'Wyckoff shuffle encountered in {prototype_label}, '
                                f'{label_wyckoff_letter} -> {equivalent_atom_set.wyckoff_letter}')
                break

    return equivalent_atom_set_list

def get_stoich_reduced_list_from_prototype(prototype_label: str) -> List[int]:
    """
    Get numerical list of stoichiometry from prototype label, i.e. "AB3_...." -> [1,3]

    Args:
        prototype_label:
            AFLOW prototype label

    Returns:
        List of reduced stoichiometric numbers
    """                        
    stoich_reduced_formula = prototype_label.split("_")[0]
    stoich_reduced_list=[]
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
            stoich_reduced_curr*=10 # will throw an error if we haven't encountered an alphabetical letter, good
            stoich_reduced_curr+=int(char)
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
    prototype_label_split = prototype_label.split('_')
    for species_wyckoff_string in prototype_label_split[3:]:
        expanded_wyckoff_letters.append('')
        curr_wyckoff_count = 0
        for char in species_wyckoff_string:
            if isalpha(char):
                if curr_wyckoff_count == 0:
                    curr_wyckoff_count = 1
                expanded_wyckoff_letters[-1] += char*curr_wyckoff_count
                curr_wyckoff_count = 0
            else:
                assert isdigit(char)
                curr_wyckoff_count *= 10 # if it's zero, we're all good
                curr_wyckoff_count += int(char)               
    return expanded_wyckoff_letters   

def get_space_group_number_from_prototype(prototype_label: str) -> int:
    return int(prototype_label.split('_')[2])

def get_pearson_symbol_from_prototype(prototype_label: str) -> str:
    return prototype_label.split('_')[1]

def get_centering_from_prototype(prototype_label: str) -> str:
    return get_pearson_symbol_from_prototype(prototype_label)[1]

def get_centering_divisor_from_prototype(prototype_label: str) -> int:
    """
    Get number of lattice points per conventional (hexagonal for rhombohedral) unit cell
    """
    return CENTERING_DIVISORS[get_centering_from_prototype(prototype_label)]

def read_shortnames() -> Dict:
    """
    This function parses ``README_PROTO.TXT``. It finds each line that (after stripping whitespace) starts with ``ANRL Label``. These are headers of sections of prototype listings. 
    It finds the column of the word ``notes``. This will be the column that the shortnames are in. 
    Skipping various non-prototype lines, the first column in each prototype line (before the ``.``) is the prototype, while the end of the line starting from the ``notes`` column, 
    cleaned up to remove whitespace and end-of-shortname comments (i.e. ``(part 3)``), is the shortname.
    
    Returns:
        A dictionary where the keys are the prototype strings, and the values are the shortnames found in the corresponding lines.
    """
    shortnames = {}
    shortname_file = "data/README_PROTO.TXT"
    notes_index = None
    with open(os.path.dirname(os.path.realpath(__file__))+'/'+shortname_file, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("ANRL Label"):
                try:
                    notes_index = line.index("notes")
                    continue
                except:
                    print("ERROR: ANRL Label line without notes header")
                    print(line)
                    sys.exit()
            # Skip this line if it's before the first ANRL label
            if notes_index == None:
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
                dum = line.split(" ")[1]
            except:
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

def get_real_to_virtual_species_map(input: Union[List[str],Atoms]) -> Dict:
    """
    Map real species to virtual species according to (alphabetized) AFLOW convention, e.g.
    for SiC return {'C':'A','Si':'B'}
    """
    if isinstance(input,Atoms):
        species = sorted(list(set(input.get_chemical_symbols())))
    else:
        species = input
    
    real_to_virtual_species_map = {}
    for i,symbol in enumerate(species):
        real_to_virtual_species_map[symbol]=chr(65+i)
    
    return real_to_virtual_species_map

def solve_for_free_params(atoms: Atoms, equation_set_list: List[EquivalentEqnSet], prototype_label: str, max_resid: float = 1e-5) -> Optional[Dict]:
    """
    Match all positions in ``atoms`` to an equation in ``equation_set_list`` to solve for the free parameters
    """
    position_set_list = group_positions_by_wyckoff(atoms,prototype_label)
    real_to_virtual_species_map = get_real_to_virtual_species_map(atoms)
    if len(position_set_list) != len(equation_set_list):
        raise inconsistentWyckoffException('Number of equivalent positions detected in Atoms object did not match the number of equivalent equations given')

    space_group_number = get_space_group_number_from_prototype(prototype_label)
    free_params_dict = {}
    position_set_matched_list = [False]*len(position_set_list)
    
    for equation_set in equation_set_list:
        # Because both equations and positions are sorted by species and wyckoff letter, this should
        # be pretty efficient
        matched_this_equation_set = False
        for i,position_set in enumerate(position_set_list):
            if position_set_matched_list[i]:
                continue
            if real_to_virtual_species_map[position_set.species] != equation_set.species:
                continue
            if not are_in_same_wyckoff_set(equation_set.wyckoff_letter,position_set.wyckoff_letter,space_group_number):
                continue
            for coeff_matrix, const_terms in zip(equation_set.coeff_matrix_list,equation_set.const_terms_list):
                for frac_position in position_set.frac_position_list:
                    possible_shifts = (-1,0,1)
                    # explore all possible shifts around zero to bring back in cell. 
                    # TODO: if this is too slow (27 possibilities), write an algorithm to determine which shifts are possible
                    for shift_list in [(x,y,z) for x in possible_shifts for y in possible_shifts for z in possible_shifts]:
                        shift_array = np.asarray(shift_list).reshape(3,1)
                        candidate_param_values,resid,_,_ = np.linalg.lstsq(coeff_matrix,frac_position-const_terms-shift_array)
                        if len(resid) == 0 or np.max(resid) < max_resid:
                            assert len(candidate_param_values) == len(equation_set.param_names)
                            for param_name,param_value in zip(equation_set.param_names,candidate_param_values):
                                assert param_name not in free_params_dict
                                free_params_dict[param_name] = param_value[0] % 1 # wrap to [0,1)
                            # should only need one to match to check off this Wyckoff position
                            position_set_matched_list[i] = True
                            matched_this_equation_set = True
                            break
                        # end loop over shifts
                    if matched_this_equation_set: break
                    # end loop over positions within a position set
                if matched_this_equation_set: break
                # end loop over equations within an equation set
            if matched_this_equation_set: break
            # end loop over position sets
        # end loop over equation sets
            
    
    if not all(position_set_matched_list):
        return None
    else:
        return free_params_dict    

class AFLOW:
    """
    Class enabling access to the AFLOW executable

    Attributes:
        aflow_executable (str): Name of the AFLOW executable
        aflow_work_dir (str): Path to the work directory
        np (int): Number of processors to use, passed to the AFLOW executable using the ``--np=...`` argument

    """

    class tooSymmetricException(Exception):
        """
        Raised when ``aflow --proto=...`` detects that the parameters requested indicate a higher symmetry
        """

    def __init__(self, aflow_executable:str="aflow", aflow_work_dir:str="",np:int=1):
        """
        Args:
            aflow_executable: Sets :attr:`aflow_executable`
            aflow_work_dir: Sets :attr:`aflow_work_dir`
            np: Sets :attr:`np`
        """
        self.aflow_executable = aflow_executable
        self.np=np
        if aflow_work_dir != "" and not aflow_work_dir.endswith("/"):
            self.aflow_work_dir = aflow_work_dir + "/"
        else:
            self.aflow_work_dir = aflow_work_dir

    def aflow_command(self, cmd: Union[str,List[str]],verbose=True) -> str:
        """
        Run AFLOW executable with specified arguments and return the output, possibly multiple times piping outputs to each other     

        Args:
            cmd: List of arguments to pass to each AFLOW executable. If it's longer than 1, multiple commands will be piped to each other            
            verbose: Whether to echo command to log file

        Raises:
            tooSymmetricException: if an ``aflow --proto=`` command complains that 
                ``the structure has a higher symmetry than indicated by the label`` 
        
        Returns:
            Output of the AFLOW command
        """
        if not isinstance(cmd,list):
            cmd = [cmd]
        
        cmd_list = [self.aflow_executable + " --np=" + str(self.np) + " " + cmd_inst
            for cmd_inst in cmd]
        cmd_str = " | ".join(cmd_list)
        if verbose:
            logger.info(cmd_str)
        try:
            return subprocess.check_output(cmd_str, shell=True, stderr=subprocess.PIPE,encoding="utf-8")
        except subprocess.CalledProcessError as exc:
            if "--proto=" in cmd_str and "The structure has a higher symmetry than indicated by the label. The correct label and parameters for this structure are:" in str(exc.stderr):
                raise self.tooSymmetricException("WARNING: the following command refused to write a POSCAR because it detected a higher symmetry: %s"%cmd_str)
            else:
                raise RuntimeError("ERROR: unexpected error from aflow command %s , error code = %d\nstderr: %s" % (cmd_str, exc.returncode, exc.stderr))
    
    def write_poscar(self, prototype_label: str, output_file: Union[str,None]=None, free_params: Union[List[float],None]=None, verbose: bool=True):
        """
        Run the ``aflow --proto`` command to write a POSCAR coordinate file corresponding to the provided AFLOW prototype designation

        Args:
            prototype_label: An AFLOW prototype label, with or without an enumeration suffix, with or without specified atomic species
            output_file: Name of the output file. If not provided, the output is written to stdout
            free_params: The free parameters of the AFLOW prototype designation. If an enumeration suffix is not included in `prototype_label` and the prototype has free parameters besides `a`, this must be provided
            verbose: Whether to echo command to log file
        """
        command = " --proto=" + prototype_label
        if free_params:
            command += " --params=" + ",".join([str(param) for param in free_params])
        if output_file is not None:
            command += " > " + self.aflow_work_dir + output_file
        return self.aflow_command([command], verbose=verbose)

    def compare_materials_dir(self, materials_subdir: str, no_scale_volume: bool=True) -> List[Dict]:
        """
        Compare a directory of materials using the aflow --compare_materials -D tool

        Args:
            materials_subdir:
                Path to the directory to compare from self.aflow_work_dir
            no_scale_volume:
                If `True`, the default behavior of allowing arbitrary scaling of structures before comparison is turned off

        Returns:        
                Attributes of representative structures, their duplicates, and groups as a whole
        """
        # TODO: For efficiency, it is possible to --add_aflow_prototype_designation to the representative structures
        # This does not help if we need duplicate prototypes (for refdata), nor for library protos (as they are not ranked by match like we need)
        command = " --compare_materials -D "
        command += self.aflow_work_dir+materials_subdir
        if no_scale_volume:
            command += " --no_scale_volume"
        output=self.aflow_command([
            command + " --screen_only --quiet --print=json"
            ])
        res_json = json.loads(output)
        return res_json

    def get_aflow_version(self)-> str:
        """
        Run the ``aflow --version`` command to get the aflow version

        Returns:
            aflow++ executable version
        """
        command = " --version"        
        output = self.aflow_command([command])
        return output.strip().split()[2]

    def compare_to_prototypes(self, input_file: str) -> List[Dict]:
        """
        Run the ``aflow --compare2prototypes`` command to compare the input structure to the AFLOW library of curated prototypes

        Args:
            input_file: path to the POSCAR file containing the structure to compare

        Returns:
            JSON list of dictionaries containing information about matching prototypes. In practice, this list should be of length zero or 1
        """

        output = self.aflow_command([            
            " --prim < " + self.aflow_work_dir + input_file,
            " --compare2prototypes --catalog=anrl --quiet --print=json"
        ])
        res_json = json.loads(output)
        return res_json
    
    def get_prototype(self,input_file: str, verbose: bool=False) -> Dict:
        """
        Run the ``aflow --prototype`` command to get the AFLOW prototype designation of the input structure

        Args:
            input_file: path to the POSCAR file containing the structure to analyze
            verbose: Whether to echo command to log file
            
        Returns:
            JSON dictionaries describing the AFLOW prototype designation (label and parameters) of the input structure.
       
        """
        output=self.aflow_command([
            " --prim < " + self.aflow_work_dir + input_file,
            " --prototype --print=json"
            ],verbose=verbose)
        res_json = json.loads(output)
        return res_json    

    def get_library_prototype_label_and_shortname(self, poscar_file: str,shortnames: Dict = read_shortnames()) -> Tuple[Union[str,None],Union[str,None]]:
        """
        Use the aflow command line tool to determine the library prototype label for a structure and look up its human-readable shortname.
        In the case of multiple results, the enumeration with the smallest misfit that is in the prototypes list is returned. If none
        of the results are in the matching prototypes list, then the prototype with the smallest misfit is returned.

        Args:
            poscar_file:
                Path to input coordinate file
            shortnames:
                Dictionary with library prototype labels as keys and human-readable "shortnames" as values.

        Returns:
            * The library prototype label for the provided compound.
            * Shortname corresponding to this prototype
        """

        comparison_results = self.compare_to_prototypes(poscar_file)
        if len(comparison_results) > 1:
            # If zero results are returned it means the prototype is not in the encyclopedia at all        
            # Not expecting a case where the number of results is greater than 1.
            raise RuntimeError(
                "{} results returned from comparison instead of zero or one as expected".format(
                    len(comparison_results)
                )
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
    
    def get_pointgroup_crystal(self, input_file: str, verbose: bool=False) -> List[Dict]:
        """
        Run the ``aflow --pointgroup_crystal`` command to get the point group operations of the provided coordinate file

        Args:
            input_file: path to the POSCAR file containing the structure to analyze

        Returns:
            JSON dictionaries describing the point group of the input structure.            
            verbose: Whether to echo command to log file
        """
        return json.loads(self.aflow_command(
            [" --pointgroup_crystal --screen_only --print=json < " + self.aflow_work_dir + input_file],
            verbose=verbose))['pgroup_xtal']
    
    def get_pointgroup_crystal_from_atoms(self, atoms: Atoms, verbose: bool=False) -> List[Dict]:
        """
        Run the ``aflow --pointgroup_crystal`` command to get the point group operations of the provided atoms object

        Args:
            atoms: Atoms object containing the structure to analyze

        Returns:
            JSON dictionaries describing the point group of the input structure.            
            verbose: Whether to echo command to log file
        """
        with NamedTemporaryFile() as f:
            atoms.write(f.name,'vasp',sort=True)
            return self.get_pointgroup_crystal(f.name,verbose=verbose)
    
    def _compare_poscars(self, poscar1: str, poscar2: str) -> Dict:
        return json.loads(self.aflow_command([' --print=JSON --compare_materials=%s,%s --screen_only --quiet'%(poscar1,poscar2)],verbose=False))
            
    def _compare_Atoms(self, atoms1: Atoms, atoms2: Atoms) -> Dict:        
        with NamedTemporaryFile() as f1, NamedTemporaryFile() as f2:
            atoms1.write(f1.name,'vasp',sort=True)
            atoms2.write(f2.name,'vasp',sort=True)
            f1.seek(0)
            f2.seek(0)
            compare = self._compare_poscars(f1.name,f2.name)
        return compare

    def get_basistransformation_rotation_originshift_from_atoms(self, atoms1: Atoms, atoms2: Atoms) -> \
        Tuple[Optional[ArrayLike],Optional[ArrayLike],Optional[ArrayLike]]:
        """
        Get operations to transform atoms2 to atoms1

        Returns:
            Tuple of arrays in the order: basis transformation, rotation, origin shift
        """
        comparison_result = self._compare_Atoms(atoms1,atoms2)
        if 'structures_duplicate' in comparison_result[0] and comparison_result[0]['structures_duplicate'] != []:
            return (
                np.asarray(comparison_result[0]['structures_duplicate'][0]['basis_transformation']),
                np.asarray(comparison_result[0]['structures_duplicate'][0]['rotation']),
                np.asarray(comparison_result[0]['structures_duplicate'][0]['origin_shift'])
            )
        else:
            return None,None,None
        
    def get_basistransformation_rotation_originshift_from_poscars(self, poscar1: str, poscar2: str) -> \
        Tuple[Optional[ArrayLike],Optional[ArrayLike],Optional[ArrayLike]]:
        """
        Get operations to transform poscar2 to poscar1

        Returns:
            Tuple of arrays in the order: basis transformation, rotation, origin shift
        """        
        comparison_result = self._compare_poscars(poscar1,poscar2)
        if 'structures_duplicate' in comparison_result[0] and comparison_result[0]['structures_duplicate'] != []:
            return (
                np.asarray(comparison_result[0]['structures_duplicate'][0]['basis_transformation']),
                np.asarray(comparison_result[0]['structures_duplicate'][0]['rotation']),
                np.asarray(comparison_result[0]['structures_duplicate'][0]['origin_shift'])
            )
        else:
            return None,None,None
    
    def build_atoms_from_prototype(
            self, species: List[str], prototype_label: str, parameter_values: List[float], primitive_cell: bool = True, verbose: bool=True, proto_file:Optional[str]=None
            ) -> Atoms:
        """
        Build an atoms object from an AFLOW prototype designation
        
        Args:
            species:
                Stoichiometric species, e.g. ``['Mo','S']`` corresponding to A and B respectively for prototype label AB2_hP6_194_c_f indicating molybdenite
            prototype_label: 
                An AFLOW prototype label, without an enumeration suffix, without specified atomic species
            parameter_values: 
                The free parameters of the AFLOW prototype designation
            primitive_cell:
                Request the primitive cell
            verbose:
                Print details in the log file
            proto_file:
                Print the output of --proto to this file

        Returns:
            Object representing unit cell of the material
        """
        assert primitive_cell, 'Can only generate primitive cells for now'
        
        with NamedTemporaryFile(mode='w+') as f, (NamedTemporaryFile(mode='w+') if proto_file is None else open(proto_file,mode='w+')) as f_with_species:            
            self.write_poscar(prototype_label,f.name,parameter_values,verbose=verbose)
            f.seek(0)
            # Add line containing species
            for i,line in enumerate(f):
                f_with_species.write(line)
                if i == 4:
                    f_with_species.write(' '.join(species)+'\n')
            f_with_species.seek(0)
            atoms = ase.io.read(f_with_species.name,format='vasp')
            
        check_number_of_atoms(atoms,prototype_label,primitive_cell)

        atoms.wrap()
        
        return atoms
    
    def get_equation_sets_from_prototype(self, prototype_label: str, parameter_values: List[float]) -> \
        List[EquivalentEqnSet]:
        """
        Get the symbolic equations for the fractional positions in the unit cell of an AFLOW prototype
        
        Args:
            prototype_label: 
                An AFLOW prototype label, without an enumeration suffix, without specified atomic species
            parameter_values: 
                The free parameters of the AFLOW prototype designation
                
        Returns:
            List of EquivalentEqnSet objects

            Each EquivalentEqnSet contains:
                - species: The species of the atoms in this set.
                - wyckoff_letter: The Wyckoff letter corresponding to this set.
                - param_names: The names of the free parameters associated with this Wyckoff position.
                - coeff_matrix_list: A list of 3 x n matrices of coefficients for the free parameters.
                - const_terms_list: A list of 3 x 1 columns of constant terms in the coordinates.
        """
        equation_poscar = self.aflow_command(f'--proto={prototype_label} --params={",".join([str(param) for param in parameter_values])} --equations_only')
        
        # get a string with one character per Wyckoff position (with possible repeated letters for positions with free params)
        wyckoff_lists = get_wyckoff_lists_from_prototype(prototype_label)
        wyckoff_joined_list = "".join(wyckoff_lists)
        
        coord_lines = equation_poscar.splitlines()[7:]
        coord_iter = iter(coord_lines)
        
        space_group_number = get_space_group_number_from_prototype(prototype_label)
        
        equation_sets = []
        
        for wyckoff_letter in wyckoff_joined_list:
            species = None # have not seen a line yet, so don't know what species it is
            param_names = None # same as above. 
            coeff_matrix_list = []
            const_terms_list = []
            # the next n positions should be equivalent corresponding to this Wyckoff position
            multiplicity_per_primitive_cell = WYCKOFF_MULTIPLICITIES[space_group_number][wyckoff_letter]/get_centering_divisor_from_prototype(prototype_label)
            assert np.isclose(multiplicity_per_primitive_cell,round(multiplicity_per_primitive_cell))
            for _ in range(round(multiplicity_per_primitive_cell)):
                line_split = next(coord_iter).split()
                if species is None:
                    species = line_split[3]
                elif line_split[3] != species:
                    raise inconsistentWyckoffException(
                        f'Encountered different species within what I thought should be the lines corresponding to Wyckoff position {wyckoff_letter}\n'
                        f'Equations obtained from prototype label {prototype_label}:\n{equation_poscar}'
                    )
                # first, get the free parameters of this line
                curr_line_free_params = set() # sympy.Symbol
                coordinate_expr_list = []
                for expression_string in line_split[:3]:
                    coordinate_expr = parse_expr(expression_string)
                    curr_line_free_params.update(coordinate_expr.free_symbols)
                    coordinate_expr_list.append(coordinate_expr)
                
                # They should all have the same number, i.e. x2,y2,z2 or x14,z14, so we can just string sort them
                curr_line_free_params = list(curr_line_free_params)
                curr_line_free_params.sort(key = lambda param: str(param))
                                
                # Each line within a Wyckoff position should have the same set of free parameters
                if param_names is None:
                    param_names = curr_line_free_params
                elif param_names != curr_line_free_params:
                    raise inconsistentWyckoffException(
                        f'Encountered different free params within what I thought should be the lines corresponding to Wyckoff position {wyckoff_letter}\n'
                        f'Equations obtained from prototype label {prototype_label}:\n{equation_poscar}'
                    )
                    
                # Transform to matrices and vectors
                a,b = linear_eq_to_matrix(coordinate_expr_list,param_names)
                
                assert a.shape == (3,len(param_names))
                assert b.shape == (3,1)
                
                coeff_matrix_list.append(matrix2numpy(a,dtype=np.float64))
                const_terms_list.append(matrix2numpy(-b,dtype=np.float64))

            # Done looping over this set of equivalent positions
            equation_sets.append(EquivalentEqnSet(
                species=species,
                wyckoff_letter=wyckoff_letter,
                param_names=[str(param_name) for param_name in param_names],
                coeff_matrix_list=coeff_matrix_list,
                const_terms_list=const_terms_list,
            ))             
        
        # do some checks
        equation_sets_iter = iter(equation_sets)
        species = None
        for species_wyckoff_list in wyckoff_lists:
            species_must_change = True
            for _ in species_wyckoff_list:
                equation_set = next(equation_sets_iter)
                if species_must_change:
                    if equation_set.species == species:
                        raise inconsistentWyckoffException(
                            'The species in the equations obtained below are inconsistent with the number and multiplicity '
                            f'of Wyckoff positions in prototype label {prototype_label}\n{equation_poscar}'
                        )
                species = equation_set.species
                species_must_change = False

        return equation_sets        
    
    def confirm_unrotated_prototype_designation(      
            self,
            reference_atoms: Atoms,
            species: List[str],
            prototype_label: str,
            parameter_values: List[float],
        ) -> bool:
        """
        Check whether the provided prototype designation recreates ``reference_atoms`` as follows:
        When the cells are in :func:`ase.cell.Cell.standard_form()`, the cells are identical.
        When both crystals are rotated to standard form (rotating the cell and keeping the fractional
        coordinates unchanged), the rotation part of the mapping the two crystals to each other found by AFLOW
        is in the point group of the reference crystal (using the generated crystal would give the same
        result). In other words, the crystals are identically oriented (but possibly translated) in reference
        to their lattice vectors, which in turn must be identical up to a rotation in reference to 
        some Cartesian coordinate system.
        
        Args:
            species:
                Stoichiometric species, e.g. ``['Mo','S']`` corresponding to A and B respectively for prototype label AB2_hP6_194_c_f indicating molybdenite
            prototype_label: 
                An AFLOW prototype label, without an enumeration suffix, without specified atomic species
            parameter_values: 
                The free parameters of the AFLOW prototype designation

        Returns:
            Whether or not the crystals match
        """
        test_atoms = self.build_atoms_from_prototype(species,prototype_label,parameter_values)
        
        if not np.allclose(reference_atoms.get_cell_lengths_and_angles(),test_atoms.get_cell_lengths_and_angles(),atol=1e-4):
            logger.info(f"Cell lengths and angles do not match.\nOriginal: {reference_atoms.get_cell_lengths_and_angles()}\n"
                        f"Regenerated: {test_atoms.get_cell_lengths_and_angles()}")
            return False
        else:
            cell_lengths_and_angles = reference_atoms.get_cell_lengths_and_angles()
        
        reference_atoms_copy = reference_atoms.copy()
        
        test_atoms.set_cell(Cell.fromcellpar(cell_lengths_and_angles),scale_atoms=True)
        reference_atoms_copy.set_cell(Cell.fromcellpar(cell_lengths_and_angles),scale_atoms=True)
        
        # the rotations below are Cartesian.
        
        _,cart_rot,_ = self.get_basistransformation_rotation_originshift_from_atoms(test_atoms,reference_atoms_copy)
        
        if cart_rot is None:
            logger.info("AFLOW failed to match the crystals")
            return False
        
        point_group_ops = self.get_pointgroup_crystal_from_atoms(reference_atoms_copy)
        
        for op in point_group_ops:
            if np.allclose(cart_rot,op['Uc'],atol=1e-4):
                logger.info("Found matching rotation")
                return True
        
        logger.info("No matching rotation found")
        return False