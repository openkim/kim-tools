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
from ase import Atoms
from ase.calculators.kim.kim import KIM
from ase.calculators.calculator import Calculator
from typing import Any, Optional, List, Union
from abc import ABC, abstractmethod
from kim_property import kim_property_create, kim_property_modify, kim_property_dump
import kim_edn
from crystal_genome_util import aflow_util
from kim_query import raw_query
from tempfile import NamedTemporaryFile
import os

__version__ = "0.1.0"
__author__ = ["ilia Nikiforov", "Eric Fuemmeler"]
__all__ = [
    "KIMTestDriverError",
    "KIMTestDriver",
    "CrystalGenomeTestDriver"
]


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
    A KIM test

    Attributes:
        model_name:
            KIM model name to use for calculations
        model:
            ASE calculator to use for calculations        
        atoms:
            List of ASE atoms objects to use as the initial configurations or to build supercells
        filename:            
            Filename to which the EDN property instance will be written
    """

    def __init__(self, model_name: Optional[str] = None, model: Optional[Calculator] = None, 
                 atoms: Optional[Union[Atoms,List[Atoms]]] = None, filename: str = "output/results.edn"):
        """
        Args:
            model_name:
                KIM extended-id of the model. Provide this or `model`
            model:
                ASE calculator to use. Provide this or `model_name`
            atoms:
                List of ASE atoms objects to use as the initial configurations or to build supercells. 
                If a single atoms object is provided, it will be converted to a single-element list
            filename:
                Path to results.edn file to be written. The default provided is the correct path to work in
                the KIM Pipeline or KIM Developer Platform
        """
        if model_name is not None:
            if model is not None:
                raise KIMTestDriverError("Please provide either a KIM model name or an ASE calculator, not both")            
            self.model_name = model_name
            self.model = KIM(model_name)
        elif model is not None:
            self.model = model
        else:
            raise KIMTestDriverError("Please provide either a KIM model name or an ASE calculator")
        
        if isinstance(atoms,List):
            self.atoms = atoms
        else:
            self.atoms = [atoms]
        
        self.property_instances = "[]"
        self.filename = filename

    @abstractmethod
    def _calculate(self, structure_index: int, **kwargs):
        """
        Abstract calculate method

        Args:
            structure_index:
                KIM tests can loop over multiple structures (i.e. crystals, molecules, etc.). This indicates which is being used for the current calculation.
                TODO: Using an index here seems un-Pythonic, any way around it?
        """
        raise NotImplementedError("Subclasses must implement the _calculate method.")

    def _write_to_file(self):
        with open(self.filename, "w") as f:
            kim_property_dump(self.property_instances, f)

    def _validate(self):
        """
        Optional physics validation of properies, to be implemented by each sublass
        """
        pass

    def __call__(self, **kwargs):
        """
        runs test and outputs results
        """
        for i,atoms in enumerate(self.atoms):
            # TODO: this seems like a very un-Pythonic way to do this, but I can't think of another way to give the _calculate
            # function a way to handle multiple initial structures except mandating that the developer always include a loop in _calculate.
            # Just passing atoms to calculate wouldn't work, what if, for example, someone has a Crystal Genome test that works directly
            # with the symmetry-reduced description?

            # still, the most common use case is an ASE calculation with Atoms, so set the calculator here
            atoms.calc = self.model
            self._calculate(i, **kwargs)

        self._validate()
        self._write_to_file()

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
        property_instances_deserialized = kim_edn.loads(self.property_instances)
        new_instance_index = len(property_instances_deserialized) + 1
        for property_instance in property_instances_deserialized:
            if property_instance["instance-id"] == new_instance_index:
                raise KIMTestDriverError("instance-id that matches the length of self.property_instances already exists.\n"
                                  "Was self.property_instances edited directly instead of using this package?")
        self.property_instances = kim_property_create(new_instance_index, property_name, self.property_instances)

    def _add_key_to_current_property_instance(self, name: str, value: Any, units: Optional[str] = None):
        """
        TODO: Add uncertainty output

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
        current_instance_index = len(kim_edn.loads(self.property_instances))
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

        self.property_instances = kim_property_modify(self.property_instances, current_instance_index, *modify_args)

################################################################################
class CrystalGenomeTestDriver(KIMTestDriver):
    """
    A Crystal Genome KIM test

    Attributes:
        stoichiometric_species: List[str]
            List of unique species in the crystal
        prototype_label: str
            AFLOW prototype label for the crystal
        parameter_names: Union[List[str],None]
            Names of free parameters of the crystal besides 'a'. May be None if the crystal is cubic with no internal DOF.
            Should have length one less than `parameter_values_angstrom`
        parameter_values_angstrom: List[float]
            List of lists of parameter values, one inner list for each equilibrium crystal structure this test will use.
            The first element in each inner list is the 'a' lattice parameter in angstrom, the rest (if present) are
            in degrees or unitless
        library_prototype_label: List[Union[str,None]]
            List of AFLOW library prototype labels, one for each equilibrium. Entries may be `None`. 
        short_name: List[Union[List[str],None]]
            Lists of human-readable short names (e.g. "FCC") for each equilibrium, if present
        cell_cauchy_stress_eV_angstrom3: List[float]
            Cauchy stress on the cell in eV/angstrom^3 (ASE units) in [xx,yy,zz,yz,xz,xy] format
        temperature_K: float
            The temperature in Kelvin
    """

    def __init__(self,
                 model_name: Optional[str] = None, 
                 model: Optional[Calculator] = None,
                 atoms: Optional[Union[List[Atoms],Atoms]] = None,
                 filename: str = "output/results.edn",
                 stoichiometric_species: Optional[List[str]] = None,
                 prototype_label: Optional[str] = None,
                 parameter_values_angstrom: Optional[Union[List[List[float]],List[float]]] = None,
                 cell_cauchy_stress_eV_angstrom3: List[float] = [0,0,0,0,0,0],
                 temperature_K: float = 0
                 ):
        """
        Args:
            model_name:
                KIM model name to use for calculations
            model:
                ASE calculator to use for calculations     
            atoms:
                List of ASE atoms objects to use as the initial configurations or to build supercells.  (NOT YET IMPLEMENTED)
                If a single atoms object is provided, it will be converted to a single-element list
            filename:
                Path to results.edn file to be written. The default provided is the correct path to work in the KIM Pipeline or KIM Developer Platform
            stoichiometric_species:
                List of unique species in the crystal
            prototype_label:
                AFLOW prototype label for the crystal
            parameter_values_angstrom:
                List of lists of AFLOW prototype parameters for the crystal.
                a (first element, always present) is in angstroms, while the other parameters 
                (present for crystals with more than 1 DOF) are in degrees or unitless. 
                If the provided list is not nested, it will be converted to a 
                If this is omitted, the parameters will be queried for
            cell_cauchy_stress_eV_angstrom3:
                Cauchy stress on the cell in eV/angstrom^3 (ASE units) in [xx,yy,zz,yz,xz,xy] format
            temperature_K:
                The temperature in Kelvin
        """        
        # Initialize model, atoms, output file        
        super().__init__(model_name,model,atoms,filename)

        self.stoichiometric_species = stoichiometric_species
        self.prototype_label = prototype_label
        self.cell_cauchy_stress_eV_angstrom3 = cell_cauchy_stress_eV_angstrom3
        self.temperature_K = temperature_K
        
        self.library_prototype_label = [None]*len(self.atoms)
        self.short_name = [None]*len(self.atoms)        

        # Only handle expected combinations of inputs, but don't raise errors for unexpected combinations, 
        # who knows what a developer might want to do?
        if (self.atoms is not None) and (self.prototype_label is None) and (parameter_values_angstrom is None):
            # placeholders for loop
            self.parameter_values_angstrom = [None]*len(self.atoms)

            for i in range(len(self.atoms)):
                self._update_aflow_designation_from_atoms(i)
                # Rebuild atoms object
                aflow = aflow_util.AFLOW()
                self.atoms[i] = aflow.build_atoms_from_prototype(self.stoichiometric_species,self.prototype_label,self.parameter_values_angstrom[i])
                # Error check->update_aflow_designation should catch it 
                self._update_aflow_designation_from_atoms(i)
        elif (stoichiometric_species is not None) and (prototype_label is not None):
            # only run this code if atoms is None, so we don't overwrite an existing atoms object
            if (parameter_values_angstrom is None) and (model_name is not None):
                self._query_aflow_designation_from_label_and_species()
            else:                
                if not isinstance(parameter_values_angstrom[0],list):
                    self.parameter_values_angstrom = [parameter_values_angstrom]
                else:
                    self.parameter_values_angstrom = parameter_values_angstrom                    
                # For now, if this constructor is called to build atoms from a fully provided AFLOW designation, don't assign library prototypes to it
                # TODO: Think about how to handle this
                self.parameter_names = ["dummy"]*(len(self.parameter_values_angstrom[0])-1) 
                # TODO: Get the list of parameter names from prototype (preferably without re-analyzing atoms object)
            aflow = aflow_util.AFLOW()
            self.atoms = []
            for parameter_values_set_angstrom in self.parameter_values_angstrom:
                self.atoms.append(aflow.build_atoms_from_prototype(stoichiometric_species,prototype_label,parameter_values_set_angstrom))
        
    def _query_aflow_designation_from_label_and_species(self):
        """
        Query for all equilibrium parameter sets for this prototype label and species in the KIM database
        """
        # TODO: Some kind of generalized query interface for all tests, this is very hand-made
        cell_cauchy_stress_Pa = [component*1.6021766e+11 for component in self.cell_cauchy_stress_eV_angstrom3]
        query_result=raw_query(
            query={
                "meta.type":"tr",
                "property-id":"tag:staff@noreply.openkim.org,2023-02-21:property/crystal-structure-npt",
                "meta.subject.extended-id":self.model_name,
                "stoichiometric-species.source-value":{
                    "$size":len(self.stoichiometric_species),
                    "$all":self.stoichiometric_species
                },
                "prototype-label.source-value":self.prototype_label,
                "cell-cauchy-stress.si-value":cell_cauchy_stress_Pa,
                "temperature.si-value":self.temperature_K
            },
            fields={
                "a.si-value":1,
                "parameter-names.source-value":1,
                "parameter-values.source-value":1,
                "library-prototype-label.source-value":1,
                "short-name.source-value":1,
                },
            database="data", limit=0, flat='on') # can't use project because parameter-values won't always exist
        if "parameter-names.source-value" in query_result[0]:
            self.parameter_names = query_result[0]["parameter-names.source-value"] # determined by prototype-label, same for all crystals
        else:
            self.parameter_names = None

        self.parameter_values_angstrom = []
        self.library_prototype_label = []
        self.short_name = []

        for parameter_set in query_result:
            self.parameter_values_angstrom.append([parameter_set["a.si-value"]*1e10])
            if "parameter-values.source-value" in parameter_set: # has params other than a
                self.parameter_values_angstrom[-1] += parameter_set["parameter-values.source-value"]
            if "library-prototype-label.source-value" in parameter_set:
                self.library_prototype_label.append(parameter_set["library-prototype-label.source-value"])
            else:
                self.library_prototype_label.append(None)
            if "short-name.source-value" in parameter_set:
                short_name = parameter_set["short-name.source-value"]
                if not isinstance(short_name,list): # Necessary because we recently changed the property definition to be a list
                    short_name = [short_name]
                self.short_name.append(short_name)
            else:
                self.short_name.append(None)            

    def _update_aflow_designation_from_atoms(self, structure_index:int, atoms:Optional[Atoms] = None):
        """
        Update the `structure_index`-th Crystal Genome crystal description fields from the corresponding self.atoms object
        or the provided atoms object        
        """
        if atoms is None:
            atoms = self.atoms[structure_index]

        aflow = aflow_util.AFLOW()
        with NamedTemporaryFile('w',delete=False) as fp: #KDP has python3.8 which is missing the convenient `delete_on_close` option
            atoms.write(fp,format='vasp')
            fp.close()
            with open(fp.name) as f:
                proto_des = aflow.get_prototype(f.name)
                libproto,short_name = aflow.get_library_prototype_label_and_shortname(f.name,aflow_util.read_shortnames())
            os.remove(fp.name)

        self.parameter_values_angstrom[structure_index] = proto_des["aflow_prototype_params_values"]
        self.library_prototype_label[structure_index] = libproto
        if short_name is None:
            self.short_name[structure_index] = None
        else:
            self.short_name[structure_index] = [short_name]

        if self.prototype_label is None:
            # we have not analyzed a single prototype yet
            assert self.stoichiometric_species is None
            self.prototype_label = proto_des["aflow_prototype_label"]
            parameter_names = proto_des["aflow_prototype_params_list"][1:]
            if len(parameter_names) > 1:
                self.parameter_names = parameter_names
            else:
                self.parameter_names = None
            self.stoichiometric_species = sorted(list(set(atoms.get_chemical_symbols())))
        else:
            if proto_des["aflow_prototype_label"] != self.prototype_label:
                raise KIMTestDriverError("It appears that the symmetry (i.e. AFLOW prototype label) is not uniform among provided "
                                  "structures, or it has changed")
            if sorted(list(set(atoms.get_chemical_symbols()))) != self.stoichiometric_species:
                raise KIMTestDriverError("It appears that the set of species is not uniform among provided "
                                  "structures, or it has changed")

    def _add_common_crystal_genome_keys_to_current_property_instance(self, structure_index: int, write_stress: bool = False, write_temp: bool = False):
        """
        Write common Crystal Genome keys -- prototype description and, optionally, stress and temperature

        Args:
            structure_index:
                Crystal Genome tests may take multiple equilibrium crystal structures of a shared prototype label and species,
                resulting in multiple property instances of the same property(s) possibly being written. This indicates which
                one is being written.
            write_stress:
                Write the `cell-cauchy-stress` key
            write_temp:
                Write the `temperature` key
        """
        self._add_key_to_current_property_instance("prototype-label",self.prototype_label)
        self._add_key_to_current_property_instance("stoichiometric-species",self.stoichiometric_species)
        self._add_key_to_current_property_instance("a",self.parameter_values_angstrom[structure_index][0],"angstrom")
        if self.parameter_names is not None:            
            self._add_key_to_current_property_instance("parameter-names",self.parameter_names)
            self._add_key_to_current_property_instance("parameter-values",self.parameter_values_angstrom[structure_index][1:])
        if self.library_prototype_label[structure_index] is not None:
            self._add_key_to_current_property_instance("library-prototype-label",self.library_prototype_label[structure_index])
        if self.short_name[structure_index] is not None:
            self._add_key_to_current_property_instance("short-name",self.short_name[structure_index])
        
        if write_stress:
            self._add_key_to_current_property_instance("cell-cauchy-stress",self.cell_cauchy_stress_eV_angstrom3,"eV/angstrom^3")
        if write_temp:
            self._add_key_to_current_property_instance("temperature",self.temperature_K,"K")

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
