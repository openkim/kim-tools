################################################################################
#
#  CDDL HEADER START
#
#  The contents of this file are subject to the terms of the Common Development
#  and Distribution License Version 1.0 (the "License").
#
#  You can obtain a copy of the license at
#  http:# www.opensource.org/licenses/CDDL-1.0.  See the License for the
#  specific language governing permissions and limitations under the License.
#
#  When distributing Covered Code, include this CDDL HEADER in each file and
#  include the License file in a prominent location with the name LICENSE.CDDL.
#  If applicable, add the following below this CDDL HEADER, with the fields
#  enclosed by brackets "[]" replaced with your own identifying information:
#
#  Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
#
#  CDDL HEADER END
#
#  Copyright (c) 2017-2019, Regents of the University of Minnesota.
#  All rights reserved.
#
#  Contributor(s):
#     Ilia Nikiforov
#     Eric Fuemmler
#
################################################################################
"""
Helper classes for KIM Test Drivers

"""
import numpy as np
from numpy.typing import ArrayLike
from ase import Atoms
from ase.calculators.kim.kim import KIM
from ase.calculators.calculator import Calculator
from typing import Any, Optional, List, Union, Dict
from abc import ABC, abstractmethod
from kim_property import kim_property_create, kim_property_modify, kim_property_dump
import kim_edn
from crystal_genome_util import aflow_util
from kim_query import raw_query
from tempfile import NamedTemporaryFile
import os
from warnings import warn
from io import StringIO

__version__ = "0.1.0"
__author__ = ["ilia Nikiforov", "Eric Fuemmeler"]
__all__ = [
    "KIMTestDriverError",
    "KIMTestDriver",
    "get_crystal_genome_designation_from_atoms",
    "verify_unchanged_symmetry",
    "CrystalGenomeTestDriver",
    "query_crystal_genome_structures",
]

# TODO: Verify uniform symmetry

################################################################################
class KIMTestDriverError(Exception):
    def __init__(self, msg):
        # Call the base class constructor with the parameters it needs
        super(KIMTestDriverError, self).__init__(msg)
        self.msg = msg

    def __str__(self):
        return self.msg

################################################################################
class KIMTestDriver(ABC):
    """
    A base class for creating KIM Test Drivers. It has attributes that are likely
    to be useful to for most KIM tests

    Attributes:
        kim_model_name: Optional[str]
            KIM model name, absent if a non-KIM ASE calculator was provided
        _calc: Calculator
            ASE calculator
        atoms: Optional[Atoms]
            ASE atoms object
        _property_instances: str
            A string containing the serialized KIM-EDN formatted property instances.
        _cached_files: Dict
            keys: filenames to be assigned to files, values: serialized strings to dump into those files. To be used for 'file' type properties
    """

    def __init__(self, model: Union[str,Calculator]):
        """
        The base class constructor only initializes the calculator and an empty property_instances

        Args:
            model:
                ASE calculator or KIM model name to use
        """
        if isinstance(model,Calculator):
            self._calc = model
            self.kim_model_name = None
        else:
            self.kim_model_name = model
            self._calc = KIM(self.kim_model_name)
        self._cached_files = {}
        self._property_instances = "[]"

    def _setup(self, atoms: Optional[Atoms] = None, **kwargs):
        """
        Set up attributes before running calculation
        """
        self.atoms = atoms

    @abstractmethod
    def _calculate(self, **kwargs):
        """
        Abstract calculate method
        """
        raise NotImplementedError("Subclasses must implement the _calculate method.")

    def write_property_instances_to_file(self,filename="output/results.edn"):
        # Write the property instances to a file at the requested path. Also dumps any cached files to the same directory
        with open(filename, "w") as f:
            kim_property_dump(self._property_instances, f)
        for cached_file in self._cached_files:        
            with open(os.path.join(os.path.dirname(filename),cached_file),"w") as f:
                f.write(self._cached_files[cached_file])
        

    def __call__(self, atoms: Optional[Atoms] = None, **kwargs):
        """
        runs and validates test
        """
        self._setup(atoms, **kwargs)
        if self.atoms is not None:
            self.atoms.calc = self._calc        
        self._calculate(**kwargs)

    def _add_property_instance(self, property_name: str):
        """
        Initialize a new property instance to self.property_instances. It will automatically get the an instance-id
        equal to the length of self.property_instances after it is added. It assumed that if you are calling this function,
        you have been only using the simplified property functions in this class and not doing any more advanced editing
        to self.property_instance using kim_property or any other methods.

        Args:
            property_name:
                The property name, e.g. "tag:staff@noreply.openkim.org,2023-02-21:property/binding-energy-crystal" or
                "binding-energy-crystal"
        """
        # DEV NOTE: I like to use the package name when using kim_edn so there's no confusion with json.loads etc.
        property_instances_deserialized = kim_edn.loads(self._property_instances)
        new_instance_index = len(property_instances_deserialized) + 1
        for property_instance in property_instances_deserialized:
            if property_instance["instance-id"] == new_instance_index:
                raise KIMTestDriverError("instance-id that matches the length of self.property_instances already exists.\n"
                                  "Was self.property_instances edited directly instead of using this package?")
        self._property_instances = kim_property_create(new_instance_index, property_name, self._property_instances)

    def _add_key_to_current_property_instance(self, name: str, value: ArrayLike, units: Optional[str] = None):
        """
        Write a key to the last element of self.property_instances. If the value is an array,
        this function will assume you want to write to the beginning of the array in every dimension.
        This function is intended to write entire keys in one go, and should not be used for modifying
        existing keys.

        WARNING! It is the developer's responsibility to make sure the array shape matches the extent
        specified in the property definition. This method uses kim_property.kim_property_modify, and
        fills the values of array keys as slices through the last dimension. If those slices are incomplete,
        kim_property automatically initializes the other elements in that slice to zero. For example,
        consider writing coordinates to a key with extent [":",3]. The correct way to write a single atom
        would be to provide [[x,y,z]]. If you accidentally provide [[x],[y],[z]], it will fill the
        field with the coordinates [[x,0,0],[y,0,0],[z,0,0]]. This will not raise an error, only exceeding
        the allowed dimesnsions of the key will do so.

        Args:
            name:
                Name of the key, e.g. "cell-cauchy-stress"
            value:
                The value of the key. The function will attempt to convert it to a NumPy array, then
                use the dimensions of the resulting array. Scalars, lists, tuples, and arrays should work.
                Data type of the elements should be str, float, or int
            units:
                The units
        """
        value_arr = np.array(value)
        value_shape = value_arr.shape
        current_instance_index = len(kim_edn.loads(self._property_instances))
        modify_args = ["key", name]
        if len(value_shape) == 0:
            modify_args += ["source-value", value]
        else:
            def recur_dimensions(prev_indices: List[int], sub_value: np.ndarray, modify_args: list):
                sub_shape = sub_value.shape
                assert len(sub_shape) != 0, "Should not have gotten to zero dimensions in the recursive function"
                if len(sub_shape) == 1:
                    # only if we have gotten to a 1-dimensional sub-array do we write stuff
                    modify_args += ["source-value", *prev_indices, "1:%d" % sub_shape[0], *sub_value]
                else:
                    for i in range(sub_shape[0]):
                        prev_indices.append(i + 1)  # convert to 1-based indices
                        recur_dimensions(prev_indices, sub_value[i], modify_args)
                        prev_indices.pop()

            prev_indices = []
            recur_dimensions(prev_indices, value_arr, modify_args)

        if units is not None:
            modify_args += ["source-unit", units]

        self._property_instances = kim_property_modify(self._property_instances, current_instance_index, *modify_args)

    def get_property_instances(self) -> Dict:
        return kim_edn.loads(self._property_instances)

################################################################################
def get_crystal_genome_designation_from_atoms(atoms: Atoms,) -> Dict:
    """
    Get crystal genome designation from an ASE atoms object.

    Returns:
        A dictionary with the following keys:
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
                List of human-readable short names (e.g. "FCC"), if present
    """
    aflow = aflow_util.AFLOW()
    cg_des = {}
    
    with NamedTemporaryFile('w',delete=False) as fp: #KDP has python3.8 which is missing the convenient `delete_on_close` option
        atoms.write(fp,format='vasp')
        fp.close()
        with open(fp.name) as f:
            proto_des = aflow.get_prototype(f.name)
            libproto,short_name = aflow.get_library_prototype_label_and_shortname(f.name,aflow_util.read_shortnames())
        os.remove(fp.name)

    cg_des["prototype_label"] = proto_des["aflow_prototype_label"]
    cg_des["stoichiometric_species"] = sorted(list(set(atoms.get_chemical_symbols())))
    parameter_names = proto_des["aflow_prototype_params_list"][1:]
    if parameter_names == []:
        cg_des["parameter_names"] = None
    else:
        cg_des["parameter_names"] = parameter_names
    cg_des["parameter_values_angstrom"] = proto_des["aflow_prototype_params_values"]
    cg_des["library_prototype_label"] = libproto
    if short_name is None:
        cg_des["short_name"] = None
    else:
        cg_des["short_name"] = [short_name]

    return cg_des

################################################################################
def verify_unchanged_symmetry(
    reference_stoichiometric_species: List[str],
    reference_prototype_label: str,
    stoichiometric_species: List[str],
    prototype_label: str,
    loose_triclinic_and_monoclinic: bool = False,
    **kwargs
    ):
    """
    Checks if stoichiometric_species and/or prototype_label have changed. Raises an error if they have

    Args:
        loose_triclinic_and_monoclinic:
            For triclinic and monoclinic space groups (1-15), the Wyckoff letters can be assigned in a non-unique way. Therefore,
            it may be useful to only check the first three parts of the prototype label: stoichiomerty, Pearson symbol and space
            group number.
    """
    checking_full_label = True
    if (int(reference_prototype_label.split("_")[2]) < 16) and loose_triclinic_and_monoclinic: # triclinic or monoclinic space group
        checking_full_label = False

    if checking_full_label:
        if reference_prototype_label != prototype_label:
            raise KIMTestDriverError("AFLOW prototype label %s differs from reference prototype label %s" % (prototype_label,reference_prototype_label))
    else:
        if reference_prototype_label.split("_")[:3] != prototype_label.split("_")[:3]:
            raise KIMTestDriverError("AFLOW prototype label %s differs from reference prototype label %s even when ignoring Wyckoff letters"%(prototype_label,reference_prototype_label))
        
    if reference_stoichiometric_species != stoichiometric_species:
        raise KIMTestDriverError("List of stoichiometric species %s does not match reference list %s" % (stoichiometric_species,reference_stoichiometric_species))
    
################################################################################
class CrystalGenomeTestDriver(KIMTestDriver):
    """
    A Crystal Genome KIM test

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
            AFLOW library prototype label, may be `None`. 
        short_name: Optional[List[str]]
            List of human-readable short names (e.g. "FCC"), if present
        cell_cauchy_stress_eV_angstrom3: List[float]
            Cauchy stress on the cell in eV/angstrom^3 (ASE units) in [xx,yy,zz,yz,xz,xy] format
        temperature_K: float
            The temperature in Kelvin
        poscar: Optional[str]
            String to be dumped as a poscar file
    """
    def _setup(self,
               atoms: Optional[Atoms] = None,
               stoichiometric_species: Optional[List[str]] = None,
               prototype_label: Optional[str] = None,
               parameter_names: Optional[List[str]] = None,
               parameter_values_angstrom: Optional[List[float]] = None,
               library_prototype_label: Optional[str] = None,
               short_name: Optional[Union[List[str],str]] = None,
               cell_cauchy_stress_eV_angstrom3: List[float] = [0,0,0,0,0,0],
               temperature_K: float = 0,
               rebuild_atoms: bool = True,
               **kwargs
               ):
        """
        Args:
            atoms:
                ASE Atoms objects to use as the initial configuration or to build supercells. 
                If this is provided, none of the arguments that are part of the Crystal Genome 
                designation should be provided, and vice versa.
            stoichiometric_species:
                List of unique species in the crystal. Required part of the Crystal Genome designation. 
            prototype_label:
                AFLOW prototype label for the crystal. Required part of the Crystal Genome designation. 
            parameter_names:
                Names of free parameters of the crystal besides 'a'. May be None if the crystal is cubic with no internal DOF.
                Should have length one less than `parameter_values_angstrom`. Part of the Crystal Genome designation.
                May be omitted for debugging, required metadata for OpenKIM Pipeline operation.
            parameter_values_angstrom:
                List of AFLOW prototype parameters for the crystal. Required part of the Crystal Genome designation. 
                a (first element, always present) is in angstroms, while the other parameters 
                (present for crystals with more than 1 DOF) are in degrees or unitless. 
            library_prototype_label: 
                AFLOW library prototype label, may be `None`. Optional part of the Crystal Genome designation.  
            short_name: 
                List of any human-readable short names (e.g. "FCC") associated with the crystal. 
                Optional part of the Crystal Genome designation.
            cell_cauchy_stress_eV_angstrom3:
                Cauchy stress on the cell in eV/angstrom^3 (ASE units) in [xx,yy,zz,yz,xz,xy] format
            temperature_K:
                The temperature in Kelvin
            rebuild_atoms:
                Normally, if you provide an Atoms object, it will be analyzed for its symmetry-reduced AFLOW description,
                and then rebuilt so that the orientation is always consistent. This can rarely cause an error due to
                rounding resulting in a higher-symmetry crystal being created during the rebuild. If you do not care
                about having your Atoms in the standard AFLOW orientation, you can turn the rebuild off.
        """ 

        super()._setup(atoms)
        self.poscar = None
        self.stoichiometric_species = stoichiometric_species        
        self.prototype_label = prototype_label
        self.parameter_names = parameter_names
        self.parameter_values_angstrom = parameter_values_angstrom
        self.library_prototype_label = library_prototype_label
        if isinstance(short_name,str):
            self.short_name = [short_name]
        else:
            self.short_name = short_name
        self.cell_cauchy_stress_eV_angstrom3 = cell_cauchy_stress_eV_angstrom3
        self.temperature_K = temperature_K

        if (
            (self.stoichiometric_species is None) !=
            (self.prototype_label is None) !=
            (self.parameter_values_angstrom is None)
        ):
            print (self._setup.__doc__)
            raise KIMTestDriverError ("\n\nYou have provided some but not all of the required parts of the Crystal Genome designation specified in the docstring above.")

        if self.atoms is not None:
            if (
                (self.stoichiometric_species is not None) or # only need to check one required part
                ((self.short_name is not None)) or
                ((self.library_prototype_label is not None)) or
                ((self.parameter_names is not None))
            ):
                print (self._setup.__doc__)
                raise KIMTestDriverError ("\n\nYou have provided an Atoms object as well as at least one part of the Crystal Genome designation specified in the docstring above.\n"
                                          "Please provide only an Atoms object or a Crystal Genome designation, not both")  
                                  
            # Updates the Crystal Genome designation class attributes according to self.atoms
            # It checks to make sure that stoichiometric_species and prototype_label have not changed,
            # But they are both None for now, so the check is skipped            
            self._update_crystal_genome_designation_from_atoms()
            if rebuild_atoms:
                # rebuild atoms for consistent orientation
                aflow = aflow_util.AFLOW()
                self.atoms = aflow.build_atoms_from_prototype(self.stoichiometric_species,self.prototype_label,self.parameter_values_angstrom)
                # Formerly there was a check here yet again to make sure symmetry hasn't changed, but I don't think it's important
        elif self.stoichiometric_species is not None: # we've already checked that if this is not None, other required parts exist as well
            # some checks and cleanup
            if (len(self.parameter_values_angstrom) > 1) and (self.parameter_names is None):
                warn("You've provided parameter values besides `a`, but no parameter names.\n"
                     "Placeholders will be inserted for debugging.")
                self.parameter_names = ["dummy"]*(len(self.parameter_values_angstrom)-1)
            aflow = aflow_util.AFLOW()
            self.atoms = aflow.build_atoms_from_prototype(self.stoichiometric_species,self.prototype_label,self.parameter_values_angstrom)
            with StringIO() as output:
                self.atoms.write(output,format='vasp')
                self.poscar=output.getvalue()                        
        else:
            warn("You've provided neither a Crystal Genome designation nor an Atoms object.\n"
                     "I won't stop you, but you better know what you're doing!")  

    def _get_crystal_genome_designation_from_atoms_and_verify_unchanged_symmetry(
            self, atoms: Optional[Atoms] = None, loose_triclinic_and_monoclinic = False
    ) -> Dict:
        """
        Get crystal genome designation from an ASE atoms object, and check if symmetry is consistent with 
        existing symmetry in this class

        Args:
            loose_triclinic_and_monoclinic:
                For triclinic and monoclinic space groups (1-15), the Wyckoff letters can be assigned in a non-unique way. Therefore,
                it may be useful to only check the first three parts of the prototype label: stoichiomerty, Pearson symbol and space
                group number. Use this if you are getting unexpected errors for monoclinic and triclinic crystals

        Returns:
            A dictionary with the following keys:
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
                    List of human-readable short names (e.g. "FCC"), if present
        """
        if atoms is None:
            atoms = self.atoms

        crystal_genome_designation = get_crystal_genome_designation_from_atoms(atoms)
        assert ((self.stoichiometric_species is None) == (self.prototype_label is None)), "self.stoichiometric_species and self.prototype_label should either both be None, or neither"
        if self.stoichiometric_species is not None:
            verify_unchanged_symmetry(
                self.stoichiometric_species,self.prototype_label,**crystal_genome_designation,loose_triclinic_and_monoclinic=loose_triclinic_and_monoclinic)
                        
        return crystal_genome_designation
    
    def _update_crystal_genome_designation_from_atoms(self, atoms: Optional[Atoms] = None, loose_triclinic_and_monoclinic = False):
        """
        Update the Crystal Genome crystal description fields from the corresponding self.atoms object
        or the provided atoms object. Additionally, cache a poscar file to write later.

        Args:
            loose_triclinic_and_monoclinic:
                For triclinic and monoclinic space groups (1-15), the Wyckoff letters can be assigned in a non-unique way. Therefore,
                it may be useful to only check the first three parts of the prototype label: stoichiomerty, Pearson symbol and space
                group number. Use this if you are getting unexpected errors for monoclinic and triclinic crystals
        """
        if atoms is None:
            atoms = self.atoms

        # get designation and check that symmetry has not changed (symmetry will not be checked if own CG designation has not been set)
        crystal_genome_designation = self._get_crystal_genome_designation_from_atoms_and_verify_unchanged_symmetry(atoms, loose_triclinic_and_monoclinic)

        with StringIO() as output:
            atoms.write(output,format='vasp')
            self.poscar=output.getvalue()

        self.stoichiometric_species = crystal_genome_designation["stoichiometric_species"]
        self.prototype_label = crystal_genome_designation["prototype_label"]
        self.parameter_names = crystal_genome_designation["parameter_names"]
        self.parameter_values_angstrom = crystal_genome_designation["parameter_values_angstrom"]
        self.library_prototype_label = crystal_genome_designation["library_prototype_label"]
        self.short_name = crystal_genome_designation["short_name"]

    def _add_common_crystal_genome_keys_to_current_property_instance(self, write_stress: bool = False, write_temp: bool = False):
        """
        Write common Crystal Genome keys -- prototype description and, optionally, stress and temperature

        Args:
            write_stress:
                Write the `cell-cauchy-stress` key
            write_temp:
                Write the `temperature` key
        """
        self._add_key_to_current_property_instance("prototype-label",self.prototype_label)
        self._add_key_to_current_property_instance("stoichiometric-species",self.stoichiometric_species)
        self._add_key_to_current_property_instance("a",self.parameter_values_angstrom[0],"angstrom")
        if self.parameter_names is not None:            
            self._add_key_to_current_property_instance("parameter-names",self.parameter_names)
            self._add_key_to_current_property_instance("parameter-values",self.parameter_values_angstrom[1:])
        if self.library_prototype_label is not None:
            self._add_key_to_current_property_instance("library-prototype-label",self.library_prototype_label)
        if self.short_name is not None:
            self._add_key_to_current_property_instance("short-name",self.short_name)        
        if write_stress:
            self._add_key_to_current_property_instance("cell-cauchy-stress",self.cell_cauchy_stress_eV_angstrom3,"eV/angstrom^3")
        if write_temp:
            self._add_key_to_current_property_instance("temperature",self.temperature_K,"K")
        if self.poscar is not None:
            current_instance_index = len(kim_edn.loads(self._property_instances))
            filename = "instance-%d.poscar"%current_instance_index
            self._cached_files[filename] = self.poscar
            self._add_key_to_current_property_instance("coordinates-file",filename)

    def _add_property_instance_and_common_crystal_genome_keys(self, property_name: str, write_stress: bool = False, write_temp: bool = False):
        """
        Initialize a new property instance to self.property_instances. It will automatically get the an instance-id
        equal to the length of self.property_instances after it is added. Then, write the common Crystal Genome
        keys to it from the attributes of this class.
        It assumed that if you are calling this function, you have been only using the simplified property functions 
        in this class and not doing any more advanced editing to self.property_instances using kim_property or any other methods. 

        Args:
            property_name:
                The property name, e.g. "tag:staff@noreply.openkim.org,2023-02-21:property/binding-energy-crystal" or
                "binding-energy-crystal"
            write_stress:
                Write the `cell-cauchy-stress` key
            write_temp:
                Write the `temperature` key
        """        
        super()._add_property_instance(property_name)
        self._add_common_crystal_genome_keys_to_current_property_instance(write_stress,write_temp)
 
################################################################################
def query_crystal_genome_structures(
            kim_model_name: str,
            stoichiometric_species: List[str],
            prototype_label: str,
            cell_cauchy_stress_eV_angstrom3: List[float] = [0,0,0,0,0,0],
            temperature_K: float = 0,
        ) -> List[Dict]:
    """
    Query for all equilibrium parameter sets for this prototype label and species in the KIM database.

    Args:
        kim_model_name: str
            KIM model name
        stoichiometric_species:
            List of unique species in the crystal. Required part of the Crystal Genome designation. 
        prototype_label:
            AFLOW prototype label for the crystal. Required part of the Crystal Genome designation. 
        cell_cauchy_stress_eV_angstrom3:
            Cauchy stress on the cell in eV/angstrom^3 (ASE units) in [xx,yy,zz,yz,xz,xy] format
        temperature_K:
            The temperature in Kelvin

    Returns:
        A list of dictionaries with the following keys:
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
                List of human-readable short names (e.g. "FCC"), if present
    """
    # TODO: Some kind of generalized query interface for all tests, this is very hand-made
    cell_cauchy_stress_Pa = [component*1.6021766e+11 for component in cell_cauchy_stress_eV_angstrom3]
    query_result=raw_query(
        query={
            "meta.type":"tr",
            "property-id":"tag:staff@noreply.openkim.org,2023-02-21:property/crystal-structure-npt",
            "meta.subject.extended-id":kim_model_name,
            "stoichiometric-species.source-value":{
                "$size":len(stoichiometric_species),
                "$all":stoichiometric_species
            },
            "prototype-label.source-value":prototype_label,
            "cell-cauchy-stress.si-value":cell_cauchy_stress_Pa,
            "temperature.si-value":temperature_K
        },
        fields={
            "a.si-value":1,
            "parameter-names.source-value":1,
            "parameter-values.source-value":1,
            "library-prototype-label.source-value":1,
            "short-name.source-value":1,
            },
        database="data", limit=0, flat='on') # can't use project because parameter-values won't always exist

    list_of_cg_des = []

    for parameter_set in query_result:
        curr_cg_des = {}        
        curr_cg_des["stoichiometric_species"] = stoichiometric_species # This was part of the query, but we provide it as output for a complete designation
        curr_cg_des["prototype_label"] = prototype_label # This was part of the query, but we provide it as output for a complete designation
        if "parameter-names.source-value" in parameter_set:
            curr_cg_des["parameter_names"] = parameter_set["parameter-names.source-value"]
        else:
            curr_cg_des["parameter_names"] = None
        # first element of parameter_values_angstrom is always present and equal to `a`
        curr_cg_des["parameter_values_angstrom"] = [parameter_set["a.si-value"]*1e10]

        if "parameter-values.source-value" in parameter_set: # has params other than a
            curr_cg_des["parameter_values_angstrom"] += parameter_set["parameter-values.source-value"]
            
        if "library-prototype-label.source-value" in parameter_set:
            curr_cg_des["library_prototype_label"] = parameter_set["library-prototype-label.source-value"]
        else:
            curr_cg_des["library_prototype_label"] = None

        if "short-name.source-value" in parameter_set:
            short_name = parameter_set["short-name.source-value"]
            if not isinstance(short_name,list): # Necessary because we recently changed the property definition to be a list
                short_name = [short_name]                
            curr_cg_des["short_name"] = short_name
        else:
            curr_cg_des["short_name"] = None
        list_of_cg_des.append(curr_cg_des)

    return list_of_cg_des

# TODO: Can probably move to different location as its not ASE related 
class Registry:
    r"""Class for registry object which acts as central source of truth."""
    __entries__ = {
        "test-drivers": {},
        "verification-checks": {},
    }

    @classmethod
    def register_test_driver(cls, name):
        def wrap(func):
            cls.__entries__["test-drivers"][name] = func
            return func

        return wrap

    @classmethod
    def register_verification_check(cls, name):
        def wrap(func):
            cls.__entries__["verification-checks"][name] = func
            return func

        return wrap

    @classmethod
    def get_test_driver_class(cls, name):
        return cls.__entries__["test-drivers"].get(name, None)

    @classmethod
    def get_verification_check_class(cls, name):
        return cls.__entries__["verification-checks"].get(name, None)

registry = Registry()
        
# If called directly, do nothing
if __name__ == "__main__":
    pass

