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
from ase.calculators.calculator import Calculator
from ase import constraints
from ase.spacegroup import symmetrize
if hasattr(constraints,'FixSymmetry'):
    from ase.constraints import FixSymmetry
elif hasattr(symmetrize,'FixSymmetry'):
    from ase.spacegroup.symmetrize import FixSymmetry
else:
    raise ImportError("Can't find `FixSymmetry` in either `ase.constraints` or `ase.spacegroup.symmetrize`")
from typing import Any, Optional, List, Union, Dict, IO
from ase.optimize import LBFGSLineSearch
from ase.optimize.optimize import Optimizer
from ase.constraints import ExpCellFilter, UnitCellFilter
from abc import ABC, abstractmethod
from kim_property import kim_property_create, kim_property_modify, kim_property_dump, get_properties, get_property_id_path
from kim_property.modify import STANDARD_KEYS_SCLAR_OR_WITH_EXTENT
import kim_edn
from .. import aflow_util
from kim_query import raw_query
from tempfile import NamedTemporaryFile
import os
from warnings import warn
from io import StringIO
import shutil
from pathlib import Path
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename='kim-tools.log',level=logging.INFO,force=True)

__author__ = ["ilia Nikiforov", "Eric Fuemmeler"]
__all__ = [
    "KIMTestDriverError",
    "KIMTestDriver",
    "get_crystal_genome_designation_from_atoms",
    "verify_unchanged_symmetry",
    "CrystalGenomeTestDriver",
    "query_crystal_genome_structures",
    "minimize_wrapper",
]    

FMAX_INITIAL = 1e-5 # Force tolerance for the optional initial relaxation of the provided cell
MAXSTEPS_INITIAL = 10000 # Maximum steps for the optional initial relaxation of the provided cell

PROP_SEARCH_PATHS_INFO=(\
'- $KIM_PROPERTY_PATH (expanding globs including recursive **)\n'
'- $PWD/local-props/**/\n'
'- $PWD/local_props/**/')

def minimize_wrapper(supercell:Atoms, fmax:float=1e-5, steps:int=10000, \
                         variable_cell:bool=True, logfile:Optional[Union[str,IO]]='-',
                         algorithm: Optimizer = LBFGSLineSearch, 
                         CellFilter: UnitCellFilter = ExpCellFilter,
                         fix_symmetry: bool = False,
                         opt_kwargs: Dict = {},
                         flt_kwargs: Dict = {}) -> None:
    """
    Use LBFGSLineSearch (default) to Minimize cell energy with respect to cell shape and
    internal atom positions.

    LBFGSLineSearch convergence behavior is as follows:
    
    - The solver returns True if it is able to converge within the optimizer
      iteration limits (which can be changed by the `steps` argument passed
      to `run`), otherwise it returns False.
    - The solver raises an exception in situations where the line search cannot
      improve the solution, typically due to an incompatibility between the
      potential's values for energy, forces, and/or stress.

    This routine attempts to minimizes the energy until the force and stress
    reduce below specified tolerances given a provided limit on the number of
    allowed steps. The code returns when convergence is achieved or no
    further progress can be made, either due to reaching the iteration step
    limit, or a stalled minimization due to line search failures.

    Args:
        supercell:
            Atomic configuration to be minimized.
        fmax:
            Force convergence tolerance (the magnitude of the force on each
            atom must be less than this for convergence)
        steps:
            Maximum number of iterations for the minimization
        variable_cell:
            True to allow relaxation with respect to cell shape
        logfile:
            Log file. `'-'` means STDOUT
        algorithm:
            ASE optimizer algorithm
        CellFilter:
            Filter to use if variable_cell is requested
        fix_symmetry:
            Whether to fix the crystallographic symmetry
        opt_kwargs:
            Dictionary of kwargs to pass to optimizer
        flt_kwargs:
            Dictionary of kwargs to pass to filter (e.g. `scalar_pressure`)
    """
    if fix_symmetry:
        symmetry = FixSymmetry(supercell)
        supercell.set_constraint(symmetry)
    if variable_cell:
        supercell_wrapped = CellFilter(supercell, **flt_kwargs)
        opt = algorithm(supercell_wrapped, logfile=logfile, **opt_kwargs)
    else:
        opt = algorithm(supercell, logfile=logfile, **opt_kwargs)
    try:
        converged = opt.run(fmax=fmax, steps=steps)
        iteration_limits_reached = not converged
        minimization_stalled = False
    except Exception as e:
        minimization_stalled = True
        iteration_limits_reached = False
        print()
        print("The following exception was caught during minimization:")
        print(repr(e))
        print()

    print("Minimization "+
        ("stalled" if minimization_stalled else "stopped" if iteration_limits_reached else "converged")+
        " after "+
        (("hitting the maximum of "+str(steps)) if iteration_limits_reached else str(opt.nsteps))+
        " steps.")
    
    if minimization_stalled or iteration_limits_reached:
        print()
        print("Final forces:")
        print(supercell.get_forces())
        print()
        print("Final stress:")
        print(supercell.get_stress())
        print()

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
            from ase.calculators.kim.kim import KIM
            self.kim_model_name = model
            self._calc = KIM(self.kim_model_name)
        self._cached_files = {}
        self._property_instances = "[]"

    def _setup(self, atoms: Optional[Atoms] = None, optimize: bool = False, **kwargs):
        """
        Set up attributes before running calculation
        """
        self.atoms = atoms
        if optimize:
            if self.atoms is None:
                raise KIMTestDriverError("You have asked to optimize the initial configuration, but did not provide an Atoms object.")        
            self.atoms.calc = self._calc        
            print("Performing minimization of initial cell...")
            print()            
            minimize_wrapper(self.atoms, fmax=FMAX_INITIAL, steps=MAXSTEPS_INITIAL, variable_cell=True)
            print()
            print("Minimized fractional positions:")
            print(self.atoms.get_scaled_positions())
            print()
            print("Minimized cell parameters:")
            print(self.atoms.cell)
            print()

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

    def __call__(self, atoms: Optional[Atoms] = None, optimize: bool = False, **kwargs):
        """
        Main operation of a Test Driver:
        
            * Run :func:`~KIMTestDriver._setup` (the base class provides a barebones version, derived classes may override)
            * If :attr:`~KIMTestDriver.atoms` is defined, set its :atttr:`~ase.atoms.Atoms.calc` to :attr:`~KIMTestDriver._calc`
            * Call :func:`~KIMTestDriver._calculate` (implemented by each individual Test Driver)
        """
        # _setup is likely overridden by 
        self._setup(atoms, optimize, **kwargs)
        if self.atoms is not None:
            self.atoms.calc = self._calc            
        # implemented by each individual Test Driver
        self._calculate(**kwargs)

    def _add_property_instance(self, property_name: str, disclaimer: Optional[str]=None):
        """
        Initialize a new property instance to self.property_instances. It will automatically get the an instance-id
        equal to the length of self.property_instances after it is added. It assumed that if you are calling this function,
        you have been only using the simplified property functions in this class and not doing any more advanced editing
        to self.property_instance using kim_property or any other methods.

        Args:
            property_name:
                The property name, e.g. "tag:staff@noreply.openkim.org,2023-02-21:property/binding-energy-crystal" or
                "binding-energy-crystal"
            disclaimer:
                An optional disclaimer commenting on the applicability of this result, e.g. 
                "This relaxation did not reach the desired tolerance."
        """
        # DEV NOTE: I like to use the package name when using kim_edn so there's no confusion with json.loads etc.
        property_instances_deserialized = kim_edn.loads(self._property_instances)
        new_instance_index = len(property_instances_deserialized) + 1
        for property_instance in property_instances_deserialized:
            if property_instance["instance-id"] == new_instance_index:
                raise KIMTestDriverError("instance-id that matches the length of self.property_instances already exists.\n"
                                  "Was self.property_instances edited directly instead of using this package?")
        existing_properties = get_properties()
        property_in_existing_properties = False
        for existing_property in existing_properties:
            if existing_property == property_name or get_property_id_path(existing_property)[3] == property_name:
                property_in_existing_properties = True

        if not property_in_existing_properties:
            print('\nThe property name or id\n%s\nwas not found in kim-properties.\n'%property_name)
            print('I will now look for an .edn file containing its definition in the following locations:\n%s\n'%PROP_SEARCH_PATHS_INFO)
            
            property_search_paths = []
            
            # environment varible
            if 'KIM_PROPERTY_PATH' in os.environ:
                property_search_paths += os.environ['KIM_PROPERTY_PATH'].split(':')
                
            # CWD
            property_search_paths.append(os.path.join(Path.cwd(),'local_props','**'))
            property_search_paths.append(os.path.join(Path.cwd(),'local-props','**'))
                    
            # recursively search for .edn files in the paths, check if they are a property definition
            # with the correct name
            
            found_custom_property = False
            
            for search_path in property_search_paths:
                if found_custom_property:
                    break
                else:
                    # hack to expand globs in both absolute and relative paths
                    if search_path[0] == '/':
                        base_path = Path('/')
                        search_glob = os.path.join(search_path[1:],'*.edn')
                    else:
                        base_path = Path()
                        search_glob = os.path.join(search_path,'*.edn')
                    
                    for path in base_path.glob(search_glob):
                        if not os.path.isfile(path): # in case there's a directory named *.edn
                            continue 
                        try:
                            path_str = str(path)
                            dict_from_edn = kim_edn.load(path_str)
                            if ('property-id') in dict_from_edn:
                                property_id = dict_from_edn['property-id']
                                if property_id == property_name or get_property_id_path(property_id)[3] == property_name:
                                    property_name = path_str
                                    found_custom_property = True
                                    break
                        except Exception as e:
                            pass
        
            if not found_custom_property:
                raise KIMTestDriverError(
                    '\nThe property name or id\n%s\nwas not found in kim-properties.\n'%property_name + \
                    'I failed to find an .edn file containing a matching "property-id" key in the following locations:\n' + PROP_SEARCH_PATHS_INFO)
        
        self._property_instances = kim_property_create(new_instance_index, property_name, self._property_instances, disclaimer)

    def _add_key_to_current_property_instance(self,
                                              name: str, 
                                              value: ArrayLike, 
                                              units: Optional[str] = None, 
                                              uncertainty_info: Optional[dict] = None):
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
            uncertainty_info:
                dictionary containing any uncertainty keys you wish to include. See https://openkim.org/doc/schema/properties-framework/
                for the possible uncertainty key names. These must be the same dimension as `value`, or they may be scalars regardless
                of the shape of `value`.
        """
        
        def recur_dimensions(prev_indices: List[int], sub_value: np.ndarray, modify_args: list, key_name: str='source-value'):
            sub_shape = sub_value.shape
            assert len(sub_shape) != 0, "Should not have gotten to zero dimensions in the recursive function"
            if len(sub_shape) == 1:
                # only if we have gotten to a 1-dimensional sub-array do we write stuff
                modify_args += [key_name, *prev_indices, "1:%d" % sub_shape[0], *sub_value]
            else:
                for i in range(sub_shape[0]):
                    prev_indices.append(i + 1)  # convert to 1-based indices
                    recur_dimensions(prev_indices, sub_value[i], modify_args, key_name)
                    prev_indices.pop()


        value_arr = np.array(value)
        value_shape = value_arr.shape

        current_instance_index = len(kim_edn.loads(self._property_instances))
        modify_args = ["key", name]
        if len(value_shape) == 0:
            modify_args += ["source-value", value]
        else:
            prev_indices = []
            recur_dimensions(prev_indices, value_arr, modify_args)

        if units is not None:
            modify_args += ["source-unit", units]

        if uncertainty_info is not None:
            for uncertainty_key in uncertainty_info:
                if not uncertainty_key in STANDARD_KEYS_SCLAR_OR_WITH_EXTENT:
                    raise KIMTestDriverError("Uncertainty key %s is not one of the allowed options %s."%(uncertainty_key,str(STANDARD_KEYS_SCLAR_OR_WITH_EXTENT)))
                uncertainty_value = uncertainty_info[uncertainty_key]
                uncertainty_value_arr = np.array(uncertainty_value)
                uncertainty_value_shape = uncertainty_value_arr.shape

                if not(len(uncertainty_value_shape) == 0 or uncertainty_value_shape == value_shape):
                    raise KIMTestDriverError("The value %s provided for uncertainty key %s has shape %s.\n"%(uncertainty_value_arr,uncertainty_key,str(uncertainty_value_shape))+\
                                             "It must either be a scalar or match the shape %s of the source value you provided."%str(value_shape))
                if len(uncertainty_value_shape) == 0:
                    modify_args += [uncertainty_key, uncertainty_value]
                else:
                    prev_indices = []
                    recur_dimensions(prev_indices, uncertainty_value_arr, modify_args, uncertainty_key)
        self._property_instances = kim_property_modify(self._property_instances, current_instance_index, *modify_args)

    def _add_file_to_current_property_instance(self,
                                              name: str, 
                                              filename: str,
                                              add_instance_index: bool = True):
        """
        add a "file" type key-value pair to the current property instance.

        Args:
            name:
                Name of the key, e.g. "restart-file"
            filename:
                The relative path to the filename. If it does not start with "output/", the file will be moved to the "output/" directory
            add_instance_index:
                By default, a numerical index will be added before the file extension or at the end of a file with no extension. This is to 
                ensure files do not get overwritten when the _calculate method is called repeatedly.
        
        Raises:
            KIMTestDriverError:
                If the provided filename does not exist
        """

        if not os.path.isfile(filename):
            raise KIMTestDriverError("Provided filename %s does not exist." % filename)
        
        if filename.split('/')[0] == 'output':
            filename_final = filename
        else:
            filename_final = os.path.join('output',filename)            

        current_instance_index = len(kim_edn.loads(self._property_instances))
        
        if add_instance_index:
            root, ext = os.path.splitext(filename_final)
            root = root + "-" + str(current_instance_index)
            filename_final = root + ext
        
        if filename_final != filename:
            shutil.move(filename,filename_final)
        
        self._property_instances = kim_property_modify(self._property_instances, current_instance_index, "key", name, "source-value", filename_final)


    @property
    def property_instances(self) -> Dict:
        return kim_edn.loads(self._property_instances)

################################################################################
def get_crystal_genome_designation_from_atoms(atoms: Atoms, get_library_prototype: bool = True, aflow_np: int = 4) -> Dict:
    """
    Get crystal genome designation from an ASE atoms object.
    
    Args:
        atoms: Atoms object to analyze
        get_library_prototype: whether to compare against prototype library
        aflow_np: Number of processors to use with AFLOW executable

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
                List of human-readable short names (e.g. "Face-Centered Cubic"), if present
    """
    aflow = aflow_util.AFLOW(np=aflow_np)
    cg_des = {}
    
    with NamedTemporaryFile('w',delete=False) as fp: #KDP has python3.8 which is missing the convenient `delete_on_close` option
        atoms.write(fp,sort=True,format='vasp')
        fp.close()
        with open(fp.name) as f:
            proto_des = aflow.get_prototype(f.name)
            (libproto,short_name) = \
                aflow.get_library_prototype_label_and_shortname(f.name,aflow_util.read_shortnames()) \
                if get_library_prototype else (None,None)
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

    Raises:
        KIMTestDriverError:
            If the symmetries of the reference and test structures are different.
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
            List of human-readable short names (e.g. "Face-Centered Cubic"), if present
        cell_cauchy_stress_eV_angstrom3: List[float]
            Cauchy stress on the cell in eV/angstrom^3 (ASE units) in [xx,yy,zz,yz,xz,xy] format
        temperature_K: float
            The temperature in Kelvin
        crystal_genome_source_structure_id: Optional[List[str]]
            Provenance identifiers of the format '[KIM test result uuid]:[instance-id]'. 
            The chains of dependencies of this test that produced this structure ends in the test listed in the test result, 
            and started with the structure computed in the specific test result and instance-id referenced. May be
            None if this test has no dependencies.
        poscar: Optional[str]
            String to be dumped as a poscar file
    """
    def _setup(self,
               atoms: Optional[Atoms] = None,
               optimize: bool = None,
               stoichiometric_species: Optional[List[str]] = None,
               prototype_label: Optional[str] = None,
               parameter_names: Optional[List[str]] = None,
               parameter_values_angstrom: Optional[List[float]] = None,
               library_prototype_label: Optional[str] = None,
               short_name: Optional[Union[List[str],str]] = None,
               cell_cauchy_stress_eV_angstrom3: List[float] = [0,0,0,0,0,0],
               temperature_K: float = 0,
               crystal_genome_source_structure_id: Optional[List[str]] = None,
               rebuild_atoms: bool = True,
               **kwargs
               ):
        """
        Args:
            atoms:
                ASE Atoms objects to use as the initial configuration or to build supercells. 
                If this is provided, none of the arguments that are part of the Crystal Genome 
                designation should be provided, and vice versa.
            optimize:
                Relax the provided Atoms object (atom positions and cell parameters). You must provide
                an Atoms object if this is True, it is not supported with a Crystal Genome designation.
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
                List of any human-readable short names (e.g. "Face-Centered Cubic") associated with the crystal. 
                Optional part of the Crystal Genome designation.
            cell_cauchy_stress_eV_angstrom3:
                Cauchy stress on the cell in eV/angstrom^3 (ASE units) in [xx,yy,zz,yz,xz,xy] format
            temperature_K:
                The temperature in Kelvin
            crystal_genome_source_structure_id:
                Provenance identifiers of the format '[KIM test result uuid]:[instance-id]'. 
                The chains of dependencies of this test that produced this structure ends in the test listed in the test result, 
                and started with the structure computed in the specific test result and instance-id referenced. May be
                None if this test has no dependencies.
            rebuild_atoms:
                Normally, if you provide an Atoms object, it will be analyzed for its symmetry-reduced AFLOW description,
                and then rebuilt so that the orientation is always consistent. This can rarely cause an error due to
                rounding resulting in a higher-symmetry crystal being created during the rebuild. If you do not care
                about having your Atoms in the standard AFLOW orientation, you can turn the rebuild off.
        """ 

        super()._setup(atoms,optimize)
        self.poscar = None
        self.stoichiometric_species = stoichiometric_species        
        self.prototype_label = prototype_label
        self.parameter_names = parameter_names
        self.parameter_values_angstrom = parameter_values_angstrom
        self.library_prototype_label = library_prototype_label
        self.crystal_genome_source_structure_id = crystal_genome_source_structure_id
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
            if optimize:
                raise KIMTestDriverError("You have asked to optimize the initial configuration while providing a Crystal Genome designation. Initial optimization is only supported when you provide an Atoms object.")            
            # some checks and cleanup
            if (len(self.parameter_values_angstrom) > 1) and (self.parameter_names is None):
                warn("You've provided parameter values besides `a`, but no parameter names.\n"
                     "Placeholders will be inserted for debugging.")
                self.parameter_names = ["dummy"]*(len(self.parameter_values_angstrom)-1)
            aflow = aflow_util.AFLOW()
            self.atoms = aflow.build_atoms_from_prototype(self.stoichiometric_species,self.prototype_label,self.parameter_values_angstrom)
            self._update_poscar()                     
        else:
            warn("You've provided neither a Crystal Genome designation nor an Atoms object.\n"
                     "I won't stop you, but you better know what you're doing!")  

    def _update_poscar(self, atoms: Optional[Atoms] = None):
        """
        Update self.poscar string from self.atoms or a provided Atoms object

        Args:
            atoms:
                The atoms object to dump, if different from ``self.atoms``
        """

        if atoms is None:
            atoms = self.atoms

        with StringIO() as output:
            atoms.write(output,format='vasp',sort=True)
            self.poscar=output.getvalue()

    def _get_crystal_genome_designation_from_atoms_and_verify_unchanged_symmetry(
            self, atoms: Optional[Atoms] = None, loose_triclinic_and_monoclinic: bool = False
    ) -> Dict:
        """
        Get Crystal Genome designation from ``self.atoms`` or a provided :class:`ase.Atoms` object, and check if symmetry is consistent with 
        existing symmetry in this class

        Args:
            atoms:
                The atoms object to analyze, if different from ``self.atoms``
            loose_triclinic_and_monoclinic:
                For triclinic and monoclinic space groups (1-15), the Wyckoff letters can be assigned in a non-unique way. Therefore,
                it may be useful to only check the first three parts of the prototype label: stoichiomerty, Pearson symbol and space
                group number. Use this if you are getting unexpected errors for monoclinic and triclinic crystals

        Returns:
            Dict:
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
                        List of human-readable short names (e.g. "Face-Centered Cubic"), if present

        Raises:
            KIMTestDriverError:
                If the symmetry of the crystal has changed

        """
        if atoms is None:
            atoms = self.atoms

        crystal_genome_designation = get_crystal_genome_designation_from_atoms(atoms)
        assert ((self.stoichiometric_species is None) == (self.prototype_label is None)), "self.stoichiometric_species and self.prototype_label should either both be None, or neither"
        if self.stoichiometric_species is not None:
            verify_unchanged_symmetry(
                self.stoichiometric_species,self.prototype_label,**crystal_genome_designation,loose_triclinic_and_monoclinic=loose_triclinic_and_monoclinic)
                        
        return crystal_genome_designation
    
    def _update_crystal_genome_designation_from_atoms(self, atoms: Optional[Atoms] = None, loose_triclinic_and_monoclinic: bool = False):
        """
        Update the Crystal Genome crystal description fields from ``self.atoms`` or a provided :class:`ase.Atoms` object. Additionally, cache a poscar file to write later.

        Args:
            atoms:
                The atoms object to analyze, if different from ``self.atoms``
            loose_triclinic_and_monoclinic:
                For triclinic and monoclinic space groups (1-15), the Wyckoff letters can be assigned in a non-unique way. Therefore,
                it may be useful to only check the first three parts of the prototype label: stoichiomerty, Pearson symbol and space
                group number. Use this if you are getting unexpected errors for monoclinic and triclinic crystals

        Raises:
            KIMTestDriverError:
                If the symmetry of the crystal has changed                
        """
        if atoms is None:
            atoms = self.atoms

        # get designation and check that symmetry has not changed (symmetry will not be checked if own CG designation has not been set)
        crystal_genome_designation = self._get_crystal_genome_designation_from_atoms_and_verify_unchanged_symmetry(atoms, loose_triclinic_and_monoclinic)

        self._update_poscar(atoms)

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
        if self.crystal_genome_source_structure_id is not None:
            self._add_key_to_current_property_instance("crystal-genome-source-structure-id",self.crystal_genome_source_structure_id)

    def _add_property_instance_and_common_crystal_genome_keys(self, property_name: str, write_stress: bool = False, write_temp: bool = False, disclaimer: Optional[str] = None):
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
            disclaimer:
                An optional disclaimer commenting on the applicability of this result, e.g. 
                "This relaxation did not reach the desired tolerance."
        """        
        super()._add_property_instance(property_name,disclaimer)
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
    This is a utility function for running the test outside of the OpenKIM pipeline. In the OpenKIM pipeline,
    this information is delivered to the test driver through the `runner` script.

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
        List[Dict]:        
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
                    List of human-readable short names (e.g. "Face-Centered Cubic"), if present
    """
    stoichiometric_species.sort()

    # TODO: Some kind of generalized query interface for all tests, this is very hand-made
    cell_cauchy_stress_Pa = [component*1.6021766e+11 for component in cell_cauchy_stress_eV_angstrom3]
    
    raw_query_args={
        "query":{
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
        "fields":{
            "a.si-value":1,
            "parameter-names.source-value":1,
            "parameter-values.source-value":1, # can't use project because parameter-values and -names won't always exist
            "library-prototype-label.source-value":1,
            "short-name.source-value":1,
            },
        "database":"data",
        "limit":0,
        "flat":"on"
    }
    
    logger.info(f"Sending below query:\n{raw_query_args}")
    
    query_result=raw_query(**raw_query_args)
    
    logger.info(f"Query result:\n{query_result}")
    
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

    print('\n!!! Found %d unique equilibrium structures from query_crystal_genome_structures() !!!\n'%len(list_of_cg_des))

    return list_of_cg_des
        
# If called directly, do nothing
if __name__ == "__main__":
    pass

