"""Tools for working with crystal prototypes using the AFLOW command line tool"""
import numpy as np
import json
import subprocess
import sys
import os
import ase
import ase.spacegroup
from ase.spacegroup.symmetrize import refine_symmetry
from curses.ascii import isalpha, isupper, isdigit
from typing import Dict, List, Tuple, Union

def get_stoich_reduced_list_from_prototype(prototype_label: str) -> List[int]:
    """
    Get numerical list of stoichiometry from prototype label, i.e. "AB3\_...." -> [1,3]

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

def get_species_list_from_string(species_string: str) -> List[str]:
    """
    Get list of chemical symbols from concatenated string of chemical symbols, i.e. "CSi" -> ["C","Si"]
    
    Args:
        species_string:
            Concatenated string of chemical symbols

    Returns:
        List of individual chemical symbols

    Raises:
        RuntimeError:
            If passed a non-alphabetical string
    """
    if any (not isalpha(character) for character in species_string):
        raise RuntimeError("Non-alphabetical character in input")

    species_list=[]
    curr_species_string=""
    for character in species_string:
        if isupper(character) and curr_species_string != "":
            species_list.append(curr_species_string)
            curr_species_string = ""
        curr_species_string+=character
    species_list.append(curr_species_string)
    return species_list

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

def get_formula_from_prototype(prototype_label: str) -> Tuple[str,int,int]:
    """
    Returns the stoichiometric formula, number of independent species in it,
    and the number of atoms per formula by analyzing the stoichiometric prefix
    in an AFLOW prototype label.

    Args:
        prototype_label :
            AFLOW prototype label

    Returns:
        * The stoichimetric formula (e.g. `AB2C3`)
        * Number of independent atoms in the stoichiometric formula (e.g. `AB2C3` has 3 independent species: `A`, `B` and `C`)
        * Number of independent atoms in the stoichiometric formula (e.g. `AB2C3` has 1 `A` + 2 `B` + 3 `C`  for 6 atoms in the formula)

    """
    # Verify that the number species matches the prototype
    formula = prototype_label.split("_")[0]
    letters = "".join([char for char in formula if not char.isdigit()])
    number_independent_species = len(letters)

    # Compute number of atoms per formula
    number_atoms_per_formula = 0
    for i in range(number_independent_species):
        i1 = formula.index(letters[i])
        if i < number_independent_species - 1:
            i2 = formula.index(letters[i + 1])
        else:
            i2 = len(formula)
        if i2 - i1 > 1:
            number_atoms_per_formula += int(formula[i1 + 1 : i2])
        else:
            number_atoms_per_formula += 1

    # Return results
    return formula, number_independent_species, number_atoms_per_formula

def get_wyckoff_info_and_cell(sgdata: Dict) -> Tuple[List[str],np.ndarray,np.ndarray]:
    """
    Parse the output from :func:`AFLOW.get_sgdata_from_prototype` to get information about Wyckoff positions and unit cell parameters
    
    Args:
        sgdata:
            JSON output of aflow --sgdata
    Returns:
        * Atomic species of the Wyckoff positions
        * Coordinates of representative Wyckoff positions
        * Conventional cell vectors
    
    """
    wyckoff_coordinates_list = []
    wyckoff_types = [] 
    for wyck in sgdata["Wyckoff_positions"]:
        wyckoff_coordinates_list.append(wyck["position"])
        wyckoff_types.append(wyck["name"])
    wyckoff_coordinates = np.array(wyckoff_coordinates_list)
    cell = sgdata["wyccar"]["lattice"]
    return wyckoff_types, wyckoff_coordinates, cell

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

    class incorrectNumAtomsException(Exception):
        """
        Raised when the number of atoms in the atoms object or in the WYCCAR returned by aflow --sgdata does not match the number of atoms in the Pearson symbol of the prototype label
        """

    class failedRefineSymmetryException(Exception):
        """
        Raised when ASE refine_symmetry function does not complete successfully
        """

    class incorrectSpaceGroupException(Exception):
        """
        Raised when spglib or aflow --sgdata detects a different space group than the one specified in the prototype label
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

    def aflow_command(self, cmd: List[str]) -> str:
        """
        Run AFLOW executable with specified arguments and return the output, possibly multiple times piping outputs to each other     

        Args:
            cmd: List of arguments to pass to each AFLOW executable. If it's longer than 1, multiple commands will be piped to each other

        Raises:
            tooSymmetricException: if an ``aflow --proto=`` command complains that 
                ``the structure has a higher symmetry than indicated by the label`` 
        
        Returns:
            Output of the AFLOW command
        """
        cmd_list = [self.aflow_executable + " --np=" + str(self.np) + cmd_inst
            for cmd_inst in cmd]
        cmd_str = " | ".join(cmd_list)                
        try:
            return subprocess.check_output(cmd_str, shell=True, stderr=subprocess.PIPE,encoding="utf-8")
        except subprocess.CalledProcessError as exc:
            if "--proto=" in cmd_str and "The structure has a higher symmetry than indicated by the label. The correct label and parameters for this structure are:" in str(exc.stderr):
                raise self.tooSymmetricException("WARNING: the following command refused to write a POSCAR because it detected a higher symmetry: %s"%cmd_str)
            else:
                raise RuntimeError("ERROR: unexpected error from aflow command %s , error code = %d\nstderr: %s" % (cmd_str, exc.returncode, exc.stderr))
    
    def write_poscar(self, prototype_label: str, output_file: Union[str,None]=None, free_params: Union[List[float],None]=None):
        """
        Run the ``aflow --proto`` command to write a POSCAR coordinate file corresponding to the provided AFLOW prototype designation

        Args:
            prototype_label: An AFLOW prototype label, with or without an enumeration suffix, with or without specified atomic species
            output_file: Name of the output file. If not provided, the output is written to stdout
            free_params: The free parameters of the AFLOW prototype designation. If an enumeration suffix is not included in `prototype_label` and the prototype has free parameters besides `a`, this must be provided
        """
        command = " --proto=" + prototype_label
        if free_params:
            command += " --params=" + ",".join([str(param) for param in free_params])
        if output_file is not None:
            command += " > " + self.aflow_work_dir + output_file
        return self.aflow_command([command])

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
    
    def get_prototype(self,input_file: str) -> Dict:
        """
        Run the ``aflow --prototype`` command to get the AFLOW prototype designation of the input structure

        Args:
            input_file: path to the POSCAR file containing the structure to analyze

        Returns:
            JSON dictionaries describing the AFLOW prototype designation (label and parameters) of the input structure.
       
        """
        output=self.aflow_command([
            " --prim < " + self.aflow_work_dir + input_file,
            " --prototype --print=json"
            ])
        res_json = json.loads(output)
        return res_json    

    def get_library_prototype_label_and_shortname(self, poscar_file: str,shortnames: Dict) -> Tuple[Union[str,None],Union[str,None]]:
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

    def get_sgdata_from_prototype(self, species: List[str], prototype_label: str, parameter_values: List[float]) -> Dict:
        """
        Without writing any files, pipe the output from aflow --prototype to aflow --sgdata to get the wyckoff info and cell

        Args:
            species:
                Stoichiometric species, e.g. ``['Mo','S']`` corresponding to A and B respectively for prototype label AB2_hP6_194_c_f indicating molybdenite
            prototype_label: 
                An AFLOW prototype label, without an enumeration suffix, without specified atomic species
            parameter_values: 
                The free parameters of the AFLOW prototype designation

        Returns:
            JSON dict containing space group information of the structure
        """
        command = [
            " --proto="+":".join([prototype_label]+species)+" --params=" + ",".join([str(param) for param in parameter_values]),
            " --sgdata --print=json"
            ]
        output = self.aflow_command(command)
        res_json = json.loads(output)
        return res_json

    def build_atoms_from_prototype(self, species: List[str], prototype_label: str, parameter_values: List[float], verbose: bool=True):
        """
        Build an atoms object from an AFLOW prototype designation
        
        Args:
            species:
                Stoichiometric species, e.g. ``['Mo','S']`` corresponding to A and B respectively for prototype label AB2_hP6_194_c_f indicating molybdenite
            prototype_label: 
                An AFLOW prototype label, without an enumeration suffix, without specified atomic species
            parameter_values: 
                The free parameters of the AFLOW prototype designation
            verbose:
                Print details

        Returns:
            Object representing conventional unit cell of the material

        Raises:
            incorrectSpaceGroupException: If space group changes during processing
            incorrectNumAtomsException: If number of atoms changes during processing
            failedRefineSymmetryException: If spglib fails

        """        
        prototype_label_list = prototype_label.split("_")
        pearson = prototype_label_list[1]
        spacegroup = int(prototype_label_list[2])

        # get the number of atoms in conventional cell from the Pearson symbol
        num_conv_cell = 0
        for character in pearson:
            if character.isdigit():
                num_conv_cell *= 10
                num_conv_cell += int(character)
        if pearson.startswith("hR"):
            num_conv_cell *= 3        

        sgdata = self.get_sgdata_from_prototype(species, prototype_label, parameter_values)
        wyckoff_types,wyckoff_coordinates,cell = get_wyckoff_info_and_cell(sgdata)
    
        if sgdata["space_group_number"]!=spacegroup:
            raise self.incorrectSpaceGroupException("WARNING: Spacegroup %d recomputed by --sgdata does not match AFLOW prototype %s"%(sgdata["space_group_number"],prototype_label))

        if len(sgdata["wyccar"]["atoms"]) != num_conv_cell:
            raise self.incorrectNumAtomsException("WARNING: Number of atoms %d recomputed by --sgdata does not match Pearson symbol of prototype %s"%(len(sgdata["wyccar"]["atoms"]),prototype_label))

        if verbose:
            print("==== Building ASE atoms object with: ====")
            print("representative atom symbols = ", wyckoff_types)
            print("representative atom coordinates = ", wyckoff_coordinates)
            print("spacegroup = ", spacegroup)
            print("cell = ", cell)
            print("=========================================")

        # Create the atoms object
        atoms = ase.spacegroup.crystal(
            symbols=wyckoff_types,
            basis=wyckoff_coordinates,
            spacegroup=spacegroup,
            cell=cell,
        )
        if len(atoms)!=num_conv_cell:
            raise self.incorrectNumAtomsException("WARNING: Number of ASE atoms %d does not match Pearson symbol of prototype %s"%(len(atoms),prototype_label))

        try:
            dataset=refine_symmetry(atoms)
        except:
            raise self.failedRefineSymmetryException("WARNING: Error while trying to refine symmetry of ASE atoms object")

        if dataset is None:
            raise self.failedRefineSymmetryException("WARNING: spglib returned None symmetry")

        if dataset["number"]!=spacegroup:
            raise self.incorrectSpaceGroupException("WARNING: spglib spacegroup %d does not match AFLOW prototype %s"%(dataset["number"],prototype_label))

        return atoms