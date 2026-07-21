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
#     Ellad B. Tadmor
#     Daniel S. Karls
#
################################################################################
"""
Helper routines for KIM Tests and Verification Checks

"""


import itertools
import logging
import math
import multiprocessing as mp
import queue as queue_module
import random
import signal
import time
import traceback
from typing import Tuple, Union

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.calculators.kim.kim import KIM
from ase.data import chemical_symbols
from ase.lattice.cubic import FaceCenteredCubic
from scipy.optimize import minimize
from scipy.signal import find_peaks

logger = logging.getLogger(__name__)
logging.basicConfig(filename="kim-tools.log", level=logging.INFO, force=True)


__all__ = [
    "KIMASEError",
    "atom_outside_cell_along_nonperiodic_dim",
    "check_if_atoms_interacting_energy",
    "check_if_atoms_interacting_force",
    "check_if_atoms_interacting",
    "get_isolated_energy_per_atom",
    "get_model_energy_cutoff",
    "get_model_species_minimum_cutoff",
    "find_equilibrium_config_FCC",
    "fractional_coords_transformation",
    "perturb_until_all_forces_sizeable",
    "randomize_positions",
    "randomize_species",
    "remove_species_not_supported_by_ASE",
    "rescale_to_get_nonzero_energy",
    "rescale_to_get_nonzero_forces",
]


################################################################################
class KIMASEError(Exception):
    def __init__(self, msg):
        # Call the base class constructor with the parameters it needs
        super(KIMASEError, self).__init__(msg)
        self.msg = msg

    def __str__(self):
        return self.msg


################################################################################
def remove_species_not_supported_by_ASE(species):
    """
    Remove any species from the 'species' list that are not supported by ASE
    """
    supported_species = chemical_symbols[1:]
    return [s for s in species if s in supported_species]


################################################################################
def randomize_species(atoms, species, seed: Union[int, None] = 13):
    """
    Given an ASE 'atoms' object, set random element for each atom selected
    from the list of available 'species' in a way that ensures all are
    represented with the same probabilities.
    """
    random.seed(seed)

    # Indefinitely iterate through species putting them at random
    # unoccupied sites until all atoms are exhausted.
    indices = list(range(len(atoms)))
    num_occupied = 0
    for element in itertools.cycle(species):
        i = random.randint(0, len(indices) - 1)
        atoms[indices[i]].symbol = element
        del indices[i]
        num_occupied += 1
        if num_occupied == len(atoms):
            break


################################################################################
def fractional_coords_transformation(cell):
    """
    Given a set of cell vectors, this will return a transformation matrix T that can be
    used to multiply any arbitrary position to get its fractional coordinates in the
    basis of those cell vectors.
    """
    simulation_cell_volume = np.linalg.det(cell)

    T = (
        np.vstack(
            (
                np.cross(cell[:, 1], cell[:, 2]),
                np.cross(cell[:, 2], cell[:, 0]),
                np.cross(cell[:, 0], cell[:, 1]),
            )
        )
        / simulation_cell_volume
    )

    return T


################################################################################
def atom_outside_cell_along_nonperiodic_dim(T, atom_coord, pbc, tol=1e-12):
    """
    Given a transformation matrix to apply to an atomic position to get its fractional
    position in the basis of some cell vectors and the corresponding boundary
    conditions, determine if the atom is outside of the cell along any non-periodic
    directions (measured using tolerance 'tol').  This is relevant when using an SM from
    a simulator such as LAMMPS that performs spatial decomposition -- atoms that leave
    the box along non-periodic dimensions are likely to become "lost."
    """
    if all(pbc):
        # Skip the actual checks
        return False

    # Calculate fractional coordinate of this atom
    atom_coord_fractional = np.dot(T, atom_coord)
    for dof in range(0, 3):
        if not pbc[dof] and (
            atom_coord_fractional[dof] < -tol or atom_coord_fractional[dof] > 1 + tol
        ):
            return True

    # If we made it to here, this atom's position is OK
    return False


################################################################################
def randomize_positions(atoms, pert_amp, seed: Union[int, None] = 13):
    """
    Given an ASE 'atoms' object, displace all atomic coordinates by a random amount in
    the range [-pert_amp, pert_amp] along *each* dimension.  Note that all atomic
    coordinates must be inside of the corresponding cell along non-periodic dimensions.
    As each atom is looped over, we continue generating perturbations until the
    displaced position is inside of the cell along any non-periodic directions
    (displacing outside the cell along periodic dimensions is allowed, although it's up
    to the calling function to wrap the positions if they need to be).
    """
    random.seed(seed)

    pbc = atoms.get_pbc()
    if all(pbc):
        # Don't need to worry about moving atoms outside of cell
        for at in range(0, len(atoms)):
            atoms[at].position += [
                pert_amp * random.uniform(-1.0, 1.0) for i in range(3)
            ]

    else:
        # Get transformation matrix to get fractional coords
        T = fractional_coords_transformation(atoms.get_cell())

        for at in range(0, len(atoms)):
            # Check if the positional coordinate of this atom is valid to begin
            # with, i.e. it's inside the box along along all non-periodic directions
            if atom_outside_cell_along_nonperiodic_dim(T, atoms[at].position, pbc):
                raise KIMASEError(
                    "ERROR: Determined that atom {} with position {} is outside of the "
                    "simulation cell ({}) along one or more non-periodic directions.  "
                    "In order to prevent atoms from being lost, they must all be "
                    "contained inside of the simulation cell along non-periodic "
                    "dimensions.".format(at, atoms[at].position, atoms.get_cell())
                )

            for dof in range(0, 3):
                coord = atoms[at].position[dof].copy()
                done = False
                while not done:
                    atoms[at].position[dof] += random.uniform(-1.0, 1.0) * pert_amp
                    if not atom_outside_cell_along_nonperiodic_dim(
                        T, atoms[at].position, pbc
                    ):
                        done = True
                    else:
                        atoms[at].position[dof] = coord


################################################################################
def get_isolated_energy_per_atom(
    model: Union[str, Calculator],
    symbol: str,
    initial_separation: float = 1.0,
    max_separation: float = 15.0,
    separation_neg_exponent: float = 4,
    quit_early_after_convergence: bool = True,
    energy_tolerance: float = 1e-12,
) -> float:
    """
    Construct a non-periodic cell containing a single atom and compute its energy.
    It tries to iteratively finetune the atomic separation for a dimer up to a
    specified precision (separation_neg_exponent). If between two successive phases
    the energy difference is less than energy_tolerance, it stops early, i.e. if
    4.0x and 4.00x are within energy_tolerance, it stops at 4.00x.
    All separations are in Angstroms.

    Args:
        model: KIM model or ASE calculator to use for calculations
        symbol: Chemical symbol
        initial_separation: Initial separation for dimer calculations
        max_separation: maximum separation for dimer calculations
        separation_neg_exponent: Number of decimal places to refine the separation
        quit_early_after_convergence: Whether to stop early if energy converges
        energy_tolerance: Energy difference tolerance for convergence check

    Returns:
        The isolated energy per atom for the requested chemical symbol
    """
    try:
        single_atom = Atoms(
            symbol,
            positions=[(0.1, 0.1, 0.1)],
            cell=(20, 20, 20),
            pbc=(False, False, False),
        )
        if isinstance(model, str):
            calc = KIM(model)
        elif isinstance(model, Calculator):
            calc = model

        single_atom.calc = calc
        energy_per_atom = single_atom.get_potential_energy()

        # Clean up KIM calculators that are able to be cleaned
        # up. This is necessary for LAMMPS ReaxFF simulators
        # encapsulated in a KIM SM,
        # as they preallocate arrays based on neighbors and
        # can crash if the lengths of neighborlists change.
        # This should only be done when a model name was passed
        # and the calculator was instantiated using KIM() in
        # this function. For example, if instead a LAMMPSLib
        # calculator was passed in directly, cleaning it
        # can break things.
        if isinstance(model, str):
            if hasattr(calc, "clean"):
                calc.clean()
            if hasattr(calc, "__del__"):
                calc.__del__()
        del single_atom

        return energy_per_atom

    except Exception:

        def _try_dimer_energy(separation):
            try:
                dimer = Atoms(
                    [symbol, symbol],
                    positions=[(0.1, 0.1, 0.1), (0.1 + separation, 0.1, 0.1)],
                    cell=(max(20, separation + 10), 20, 20),
                    pbc=(False, False, False),
                )
                if isinstance(model, str):
                    calc = KIM(model)
                elif isinstance(model, Calculator):
                    calc = model
                dimer.calc = calc

                total_energy = dimer.get_potential_energy()
                energy_per_atom = total_energy / 2.0

                if isinstance(model, str):
                    if hasattr(calc, "clean"):
                        calc.clean()
                    if hasattr(calc, "__del__"):
                        calc.__del__()
                del dimer

                return energy_per_atom

            except Exception:
                try:
                    if isinstance(model, str):
                        if hasattr(calc, "clean"):
                            calc.clean()
                        if hasattr(calc, "__del__"):
                            calc.__del__()
                    if "dimer" in locals():
                        del dimer
                except Exception:
                    pass
                return None

        # Start with integer separations: 1.0, 2.0, 3.0, 4.0, ...
        last_successful_separation = None
        last_successful_energy = None

        separation = initial_separation
        while separation <= max_separation:
            energy = _try_dimer_energy(separation)
            if energy is not None:
                last_successful_separation = separation
                last_successful_energy = energy
                separation += 1.0
            else:
                break

        if last_successful_separation is None:
            raise RuntimeError(
                f"Failed to obtain isolated energy for {symbol} - no separations worked"
            )

        # refine
        current_separation = last_successful_separation
        previous_phase_energy = last_successful_energy

        for decimal_place in range(1, separation_neg_exponent + 1):

            step_size = 10 ** (-decimal_place)

            for i in range(1, 10):
                test_sep = current_separation + step_size
                energy = _try_dimer_energy(test_sep)

                if energy is not None:
                    current_separation = test_sep
                    last_successful_energy = energy
                else:
                    break

            if quit_early_after_convergence:
                energy_diff = abs(last_successful_energy - previous_phase_energy)
                if energy_diff <= energy_tolerance:
                    return last_successful_energy

            previous_phase_energy = last_successful_energy

        return last_successful_energy


################################################################################
def rescale_to_get_nonzero_energy(atoms, isolated_energy_per_atom, etol):
    """
    If the given configuration has a potential energy, relative to the sum of the
    isolated energy corresponding to each atom present, smaller in magnitude than 'etol'
    (presumably because the distance between atoms is too large), rescale it making it
    smaller.  The 'isolated_energy_per_atom' arg should be a dict containing an entry
    for each atomic species present in the atoms object
    (additional entries are ignored).
    """
    num_atoms = len(atoms)
    if num_atoms < 2:
        # If we're only using a single atom, then we need to make sure that the cell is
        # periodic along at least one direction
        if num_atoms == 1:
            pbc = atoms.get_pbc()
            if not any(pbc):
                raise RuntimeError(
                    "ERROR: If only a single atom is present, the cell must "
                    "be periodic along at least one direction."
                )
        else:
            raise RuntimeError(
                "ERROR: Invalid configuration. Must have at least one atom"
            )

    if not isinstance(isolated_energy_per_atom, dict):
        raise ValueError(
            "Argument 'isolated_energy_per_atom' passed to "
            "rescale_to_get_nonzero_energy must be a dict containing and entry "
            "for each atomic species present in the atoms object."
        )

    # Check for any flat directions in the initial configuration and ignore them when
    # determining whether to stop rescaling
    pmin = atoms.get_positions().min(axis=0)  # minimum x,y,z coordinates
    pmax = atoms.get_positions().max(axis=0)  # maximum x,y,z coordinates
    delp = pmax - pmin  # system extent across x, y, z
    flat = [(extent <= np.finfo(extent).tiny) for extent in delp]

    # Compute the "trivial energy", i.e. the energy assuming none of the atoms interact
    # with each other at all
    species_of_each_atom = atoms.get_chemical_symbols()
    energy_trivial = 0.0
    for atom_species in species_of_each_atom:
        energy_trivial += isolated_energy_per_atom[atom_species]

    # Rescale cell and atoms
    cell = atoms.get_cell()
    energy = atoms.get_potential_energy()
    adjusted_energy = energy - energy_trivial
    if abs(adjusted_energy) <= etol:
        pmin = atoms.get_positions().min(axis=0)  # minimum x,y,z coordinates
        pmax = atoms.get_positions().max(axis=0)  # maximum x,y,z coordinates
        extent_along_nonflat_directions = [
            extent for direction, extent in enumerate(delp) if not flat[direction]
        ]
        delpmin = min(extent_along_nonflat_directions)
        while delpmin > np.finfo(delpmin).tiny:
            atoms.positions *= 0.5  # make configuration half the size
            cell *= 0.5  # make cell half the size
            atoms.set_cell(cell)  # need to adjust cell in case it's periodic
            delpmin *= 0.5
            energy = atoms.get_potential_energy()
            adjusted_energy = energy - energy_trivial
            if abs(adjusted_energy) > etol:
                return  # success!

        # Get species and write out error
        raise RuntimeError(
            "ERROR: Unable to scale configuration down to nonzero energy.  This was "
            "determined by computing the total potential energy relative to the sum of "
            "the isolated energy corresponding to each atom present and checking if "
            "the magnitude of the difference was larger than the supplied tolerance of "
            "{} eV.  This may mean that the species present in the cell ({}) do not "
            "have a non-trivial energy interaction for the  potential being used."
            "".format(etol, set(species_of_each_atom))
        )


################################################################################
def check_if_atoms_interacting_energy(model, symbols, etol):
    """
    First, get the energy of a single isolated atom of each species given in 'symbols'.
    Then, construct a dimer consisting of these two species and try to decrease its bond
    length until a discernible difference in the energy (from the sum of the isolated
    energy of each species) is detected.  The 'symbols' arg should be a list or tuple of
    length 2 indicating which species pair to check, e.g. to check if Al interacts with
    Al, one should specify ['Al', 'Al'].
    """
    if not isinstance(symbols, (list, tuple)) or len(symbols) != 2:
        raise ValueError(
            "Argument 'symbols' passed to check_if_atoms_interacting_energy "
            "must be a list of tuple of length 2 indicating the species pair to "
            "check"
        )

    isolated_energy_per_atom = {}
    isolated_energy_per_atom[symbols[0]] = get_isolated_energy_per_atom(
        model, symbols[0]
    )
    isolated_energy_per_atom[symbols[1]] = get_isolated_energy_per_atom(
        model, symbols[1]
    )

    dimer = Atoms(
        symbols,
        positions=[(0.1, 0.1, 0.1), (5.1, 0.1, 0.1)],
        cell=(20, 20, 20),
        pbc=(False, False, False),
    )
    calc = KIM(model)
    dimer.calc = calc
    try:
        rescale_to_get_nonzero_energy(dimer, isolated_energy_per_atom, etol)
        atoms_interacting = True
        return atoms_interacting
    except:  # noqa: E722
        atoms_interacting = False
        return atoms_interacting
    finally:
        if hasattr(calc, "clean"):
            calc.clean()
        if hasattr(calc, "__del__"):
            calc.__del__()
        del dimer


################################################################################
def check_if_atoms_interacting_force(model, symbols, ftol):
    """
    Construct a dimer and try to decrease its bond length until the force acting on each
    atom is larger than 'ftol' in magnitude.  The 'symbols' arg should be a list or
    tuple of length 2 indicating which species pair to check, e.g. to check if Al
    interacts with Al, one should specify ['Al', 'Al'].
    """
    if not isinstance(symbols, (list, tuple)) or len(symbols) != 2:
        raise ValueError(
            "Argument 'symbols' passed to check_if_atoms_interacting_force "
            "must be a list of tuple of length 2 indicating the species pair to "
            "check"
        )

    dimer = Atoms(
        symbols,
        positions=[(0.1, 0.1, 0.1), (5.1, 0.1, 0.1)],
        cell=(20, 20, 20),
        pbc=(False, False, False),
    )
    calc = KIM(model)
    dimer.calc = calc
    try:
        rescale_to_get_nonzero_forces(dimer, ftol)
        atoms_interacting = True
        return atoms_interacting
    except:  # noqa: E722
        atoms_interacting = False
        return atoms_interacting
    finally:
        if hasattr(calc, "clean"):
            calc.clean()
        if hasattr(calc, "__del__"):
            calc.__del__()
        del dimer


################################################################################
def check_if_atoms_interacting(
    model, symbols, check_energy=True, etol=1e-6, check_force=True, ftol=1e-3
):
    """
    Check to see whether non-trivial energy and/or forces can be detected using the
    current model.  The 'symbols' arg should be a list or tuple of length 2 indicating
    which species pair to check, e.g. to check if Al interacts with Al, one should
    specify ['Al', 'Al'].
    """
    if check_energy and not check_force:
        return check_if_atoms_interacting_energy(model, symbols, etol)
    elif not check_energy and check_force:
        return check_if_atoms_interacting_force(model, symbols, ftol)

    elif check_energy and check_force:
        atoms_interacting_energy = check_if_atoms_interacting_energy(
            model, symbols, etol
        )
        atoms_interacting_force = check_if_atoms_interacting_force(model, symbols, ftol)
        return atoms_interacting_energy, atoms_interacting_force


################################################################################
def rescale_to_get_nonzero_forces(atoms, ftol):
    """
    If the given configuration has force components which are all smaller in absolute
    value than 'ftol' (presumably because the distance between atoms is too large),
    rescale it to be smaller until the largest force component in absolute value is
    greater than or equal to 'ftol'.  In a perfect crystal, the crystal is rescaled
    until the atoms on the surface reach the minimum value (internal atoms padded with
    another atoms around them will have zero force).  Note that any periodicity is
    turned off for the rescaling and then restored at the end.
    """
    if len(atoms) < 2:
        raise KIMASEError(
            "ERROR: Invalid configuration. Must have at least 2 atoms. Number of atoms "
            "= {}".format(len(atoms))
        )

    # Check for any flat directions in the initial configuration and ignore them when
    # determining whether to stop rescaling
    pmin = atoms.get_positions().min(axis=0)  # minimum x,y,z coordinates
    pmax = atoms.get_positions().max(axis=0)  # maximum x,y,z coordinates
    delp = pmax - pmin  # system extent across x, y, z
    flat = [(extent <= np.finfo(extent).tiny) for extent in delp]

    # Temporarily turn off any periodicity
    pbc_save = atoms.get_pbc()
    cell = atoms.get_cell()
    atoms.set_pbc([False, False, False])
    # Rescale cell and atoms
    forces = atoms.get_forces()
    if np.isnan(forces).any():
        raise RuntimeError("ERROR: Computed forces include at least one nan.")
    fmax = max(abs(forces.min()), abs(forces.max()))  # find max in abs value
    if fmax < ftol:
        pmin = atoms.get_positions().min(axis=0)  # minimum x,y,z coordinates
        pmax = atoms.get_positions().max(axis=0)  # maximum x,y,z coordinates
        extent_along_nonflat_directions = [
            extent for direction, extent in enumerate(delp) if not flat[direction]
        ]
        delpmin = min(extent_along_nonflat_directions)
        while delpmin > np.finfo(delpmin).tiny:
            atoms.positions *= 0.75  # make configuration 3/4 the size
            cell *= 0.75  # make cell 3/4 the size
            delpmin *= 0.75
            forces = atoms.get_forces()  # get max force
            fmax = max(abs(forces.min()), abs(forces.max()))
            if fmax >= ftol:
                # Restore periodicity
                atoms.set_pbc(pbc_save)
                atoms.set_cell(cell)
                return  # success!
        raise KIMASEError(
            "ERROR: Unable to scale configuration down to nonzero forces."
        )
    else:
        # Restore periodicity
        atoms.set_pbc(pbc_save)


################################################################################
def perturb_until_all_forces_sizeable(
    atoms, pert_amp, minfact=0.1, maxfact=5.0, max_iter=1000
):
    """
    Keep perturbing atoms in the ASE 'atoms' object until all force components on each
    atom have an absolute value of least 'minfact' times the largest (in absolute value)
    component across all force vectors coming in.  Note that all atomic coordinates must
    be inside of the corresponding cell along non-periodic dimensions.  Perturbations
    leading to a force component on any atom that is larger than 'maxfact' times the
    largest force component coming in are rejected.  Perturbations leading to atoms
    outside of the span across x, y, and z of the atomic positions coming in or outside
    of the simulation cell along non-periodic directions are also rejected.  The process
    repeats until max_iter iterations have been reached, at which point an exception is
    raised.
    """
    pbc = atoms.get_pbc()

    # Get transformation matrix to get fractional coords
    T = fractional_coords_transformation(atoms.get_cell())

    # First, ensure that all atomic positions are inside of the cell along non-periodic
    # dimensions
    if not all(pbc):
        for at in range(0, len(atoms)):
            # Check if the positional coordinate of this atom is valid to begin
            # with, i.e. it's inside the box along along all non-periodic directions
            if atom_outside_cell_along_nonperiodic_dim(T, atoms[at].position, pbc):
                raise KIMASEError(
                    "ERROR: Determined that atom {} with position {} is outside of the "
                    "simulation cell ({}) along one or more non-periodic directions.  "
                    "In order to prevent atoms from being lost, they must all be "
                    "contained inside of the simulation cell along non-periodic "
                    "dimensions.".format(at, atoms[at].position, atoms.get_cell())
                )

    forces = atoms.get_forces()
    if np.isnan(forces).any():
        raise RuntimeError("ERROR: Computed forces include at least one nan.")
    fmax = max(abs(forces.min()), abs(forces.max()))  # find max in abs value

    if fmax <= 1e2 * np.finfo(float).eps:
        raise KIMASEError(
            "ERROR: Largest force component on configuration is "
            "less than or equal to 1e2*machine epsilon. Cannot proceed."
        )

    pmin = atoms.get_positions().min(axis=0)  # minimum x,y,z coordinates
    pmax = atoms.get_positions().max(axis=0)  # maximum x,y,z coordinates
    saved_posns = atoms.get_positions().copy()
    saved_forces = atoms.get_forces().copy()
    some_forces_too_small = True

    # Counter to enforce max iterations
    iters = 0

    while some_forces_too_small:
        for at in range(0, len(atoms)):
            for dof in range(0, 3):
                if abs(forces[at, dof]) < minfact * fmax:
                    done = False
                    coord = atoms[at].position[dof].copy()
                    while not done:
                        atoms[at].position[dof] += random.uniform(-1.0, 1.0) * pert_amp
                        if (
                            pmin[dof] <= atoms[at].position[dof] <= pmax[dof]
                        ) and not atom_outside_cell_along_nonperiodic_dim(
                            T, atoms[at].position, pbc
                        ):
                            done = True
                        else:
                            atoms[at].position[dof] = coord
        try:
            forces = atoms.get_forces()
            if np.isnan(forces).any():
                raise RuntimeError("ERROR: Computed forces include at least one nan.")
            fmax_new = max(abs(forces.min()), abs(forces.max()))
            if fmax_new > maxfact * fmax:
                # forces too large, abort perturbation
                atoms.set_positions(saved_posns)
                forces = saved_forces.copy()
                continue
            fmin_new = min(abs(forces.min()), abs(forces.max()))
            if fmin_new > minfact * fmax:
                some_forces_too_small = False
        except:  # noqa: E722
            # force calculation failed, abort perturbation
            atoms.set_positions(saved_posns)
            continue
        finally:
            iters = iters + 1
            if iters == max_iter:
                raise KIMASEError(
                    "Maximum iterations ({}) exceeded in call to "
                    "function perturb_until_all_forces_sizeable()".format(max_iter)
                )


################################################################################
def get_model_species_minimum_cutoff(
    model: Union[str, Calculator],
    species: list,
    xtol: float = 1e-8,
    etol_coarse: float = 1e-6,
    etol_fine: float = 1e-15,
    max_bisect_iters: int = 1000,
    max_upper_cutoff_bracket: float = 20.0,
) -> float:
    """
    Given a model and a list of species, construct all permutations of dimers and
    compute their energy-cutoff i.e the distance at which their interaction becomes
    non-trivial. Return the smallest such distance among all pairs of species.
    This method calls get_model_energy_cutoff
    """
    min_cutoff = np.inf
    for i in range(len(species)):
        si: str = species[i]
        Rii: float = get_model_energy_cutoff(
            model,
            [si, si],
            xtol,
            etol_coarse,
            etol_fine,
            max_bisect_iters,
            max_upper_cutoff_bracket,
        )
        min_cutoff = min(min_cutoff, Rii)
        for j in range(i + 1, len(species)):
            sj: str = species[j]
            try:
                Rij: float = get_model_energy_cutoff(
                    model,
                    [sj, si],
                    xtol,
                    etol_coarse,
                    etol_fine,
                    max_bisect_iters,
                    max_upper_cutoff_bracket,
                )
                min_cutoff = min(min_cutoff, Rij)
            except Exception:
                pass
    return min_cutoff


################################################################################
def get_model_energy_cutoff(
    model: Union[str, Calculator],
    symbols: list,
    xtol: float = 1e-8,
    etol_coarse: float = 1e-6,
    etol_fine: float = 1e-15,
    max_bisect_iters: int = 1000,
    max_upper_cutoff_bracket: float = 20.0,
) -> float:
    """
    Compute the distance at which energy interactions become non-trival for a given
    model and a species pair it supports.  This is done by constructing a dimer composed
    of these species in a large finite box, increasing the separation if necessary until
    the total potential energy is within 'etol_fine' of the sum of the corresponding
    isolated energies, and then shrinking the separation until the energy differs from
    that value by more than 'etol_coarse'.  Using these two separations to bound the
    search range, bisection is used to refine in order to locate the cutoff.  The
    'symbols' arg should be a list or tuple of length 2 indicating which species pair to
    check, e.g. to get the energy cutoff of Al with Al, one should specify ['Al', 'Al'].

    This function is based on the content of the DimerContinuityC1__VC_303890932454_002
    Verification Check in OpenKIM [1-3].

    [1] Tadmor E. Verification Check of Dimer C1 Continuity v002. OpenKIM; 2018.
        doi:10.25950/43d2c6d5

    [2] Tadmor EB, Elliott RS, Sethna JP, Miller RE, Becker CA. The potential of
        atomistic simulations and the Knowledgebase of Interatomic Models. JOM.
        2011;63(7):17. doi:10.1007/s11837-011-0102-6

    [3] Elliott RS, Tadmor EB. Knowledgebase of Interatomic Models (KIM) Application
        Programming Interface (API). OpenKIM; 2011. doi:10.25950/ff8f563a
    """
    from scipy.optimize import bisect

    def get_dimer_positions(a, large_cell_len):
        """
        Generate positions for a dimer of length 'a' centered in a finite simulation box
        with side length 'large_cell_len'
        """
        half_cell = 0.5 * large_cell_len
        positions = [
            [half_cell - 0.5 * a, half_cell, half_cell],
            [half_cell + 0.5 * a, half_cell, half_cell],
        ]
        return positions

    def energy(a, dimer, large_cell_len, einf):
        dimer.set_positions(get_dimer_positions(a, large_cell_len))
        return dimer.get_potential_energy() - einf

    def energy_cheat(a, dimer, large_cell_len, offset, einf):
        dimer.set_positions(get_dimer_positions(a, large_cell_len))
        return (dimer.get_potential_energy() - einf) + offset

    if not isinstance(symbols, (list, tuple)) or len(symbols) != 2:
        raise ValueError(
            "Argument 'symbols' passed to check_if_atoms_interacting_energy "
            "must be a list of tuple of length 2 indicating the species pair to "
            "check"
        )

    isolated_energy_per_atom = {}
    isolated_energy_per_atom[symbols[0]] = get_isolated_energy_per_atom(
        model, symbols[0]
    )
    isolated_energy_per_atom[symbols[1]] = get_isolated_energy_per_atom(
        model, symbols[1]
    )
    einf = isolated_energy_per_atom[symbols[0]] + isolated_energy_per_atom[symbols[1]]

    # First, establish the upper bracket cutoff by starting at 'b_init' Angstroms and
    # incrementing by 'db' until
    b_init = 4.0

    # Create finite box of size large_cell_len
    large_cell_len = 50
    dimer = Atoms(
        symbols,
        positions=get_dimer_positions(b_init, large_cell_len),
        cell=(large_cell_len, large_cell_len, large_cell_len),
        pbc=(False, False, False),
    )
    calc = KIM(model)
    dimer.calc = calc

    db = 2.0
    b = b_init - db
    still_interacting = True
    while still_interacting:
        b += db
        if b > max_upper_cutoff_bracket:
            if hasattr(calc, "clean"):
                calc.clean()
            if hasattr(calc, "__del__"):
                calc.__del__()

            raise KIMASEError(
                "Exceeded limit on upper bracket when determining cutoff "
                "search range"
            )
        else:
            eb = energy(b, dimer, large_cell_len, einf)
            if abs(eb) < etol_fine:
                still_interacting = False

    a = b
    da = 0.01
    not_interacting = True
    while not_interacting:
        a -= da
        if a < 0:
            if hasattr(calc, "clean"):
                calc.clean()
            if hasattr(calc, "__del__"):
                calc.__del__()

            raise RuntimeError(
                "Failed to determine lower bracket for cutoff search using etol_coarse "
                "= {}.  This may mean that the species pair provided ({}) does not "
                "have a non-trivial energy interaction for the potential being "
                "used.".format(etol_coarse, symbols)
            )
        else:
            ea = energy(a, dimer, large_cell_len, einf)
            if abs(ea) > etol_coarse:
                not_interacting = False

    # NOTE: Some Simulator Models have a history dependence due to them maintaining
    #       charges from the previous energy evaluation to use as an initial guess
    #       for the next charge equilibration.  We therefore have to treat them not
    #       as single-valued functions but as distributions, i.e.  for a given
    #       configuration you might get any of a range of energy values depending on
    #       the history of your previous energy evaluations.  This is particularly
    #       problematic for this step, where we set up a bisection problem in order
    #       to determine the cutoff radius of the model.  Our solution for this
    #       specific case is to make a very crude estimate of the variance of that
    #       distribution with a 10% factor of safety on it.
    eb_new = energy(b, dimer, large_cell_len, einf)
    eb_error = abs(eb_new - eb)

    # compute offset to ensure that energy before and after cutoff have
    # different signs
    if ea < eb:
        offset = -eb + 1.1 * eb_error + np.finfo(float).eps
    else:
        offset = -eb - 1.1 * eb_error - np.finfo(float).eps

    rcut, results = bisect(
        energy_cheat,
        a,
        b,
        args=(dimer, large_cell_len, offset, einf),
        full_output=True,
        xtol=xtol,
        maxiter=max_bisect_iters,
    )

    # General clean-up
    if hasattr(calc, "clean"):
        calc.clean()
    if hasattr(calc, "__del__"):
        calc.__del__()

    if not results.converged:
        raise RuntimeError(
            "Bisection search to find cutoff distance did not converge "
            "within {} iterations with xtol = {}".format(max_bisect_iters, xtol)
        )
    else:
        return rcut


################################################################################
# FIND-EQUILIBRIUM-CONFIGURATION-FCC
################################################################################

DEFAULT_FCC_NCELLS_PER_SIDE = 1


def _species_label(species_list: list[str]) -> str:
    return "-".join(species_list)


def _round_alat(alat: float) -> float:
    return round(float(alat), 8)


def fcc_atoms_in_supercell(ncells_per_side: int) -> int:
    """
    Calculate the total number of atoms in an FCC supercell.

    Args:
        ncells_per_side: The number of primitive unit cells aligned along
          a single axis of the supercell box.

    Returns:
        The total count of atom sites within the specified supercell as an
        integer.
    """
    return int(4 * int(ncells_per_side) ** 3)


def make_fcc_template(
    species_list: list[str],
    ncells_per_side: int = DEFAULT_FCC_NCELLS_PER_SIDE,
):
    """
    Create a generic FCC template large enough to contain at least
    len(species_list) atoms. The actual species are assigned later.

    Args:
        species_list: Chemical symbols of the target elements to be added later,
          used here strictly to determine the minimum required atom count.
        ncells_per_side: Initial number of unit cells along each axis to define
          the supercell size. Dynamically increments if the resulting supercell
          is too small. Defaults to DEFAULT_FCC_NCELLS_PER_SIDE.

    Returns:
        A tuple containing:
            - atoms: An ASE Atoms object representing the generated hydrogen-based
              Face-Centered Cubic supercell template.
            - ncells_per_side: The final integer number of unit cells per side
              used to satisfy the atom count constraint.
    """

    while True:
        atoms = FaceCenteredCubic(
            size=(ncells_per_side, ncells_per_side, ncells_per_side),
            latticeconstant=1.0,
            symbol="H",
            pbc=True,
        )

        if len(atoms) < len(species_list):
            ncells_per_side += 1
        else:
            break

    return atoms, ncells_per_side


def make_fcc_reference_config(
    species_list: list[str],
    ncells_per_side: int = DEFAULT_FCC_NCELLS_PER_SIDE,
    seed: int = 13,
) -> dict:
    """
    Create one fixed FCC reference configuration.

    All later energy evaluations use this same atom ordering, same species
    assignment, same scaled positions, and same number of unit cells. Only
    the lattice constant changes.

    Args:
        species_list: Chemical symbols of the elements to assign to the lattice
          sites.
        ncells_per_side: Number of unit cells along each axis to define the
          supercell dimensions. Defaults to DEFAULT_FCC_NCELLS_PER_SIDE.
        seed: Random seed used to determine reproducible chemical species
          ordering and fractional atom site placements. Defaults to 13.

    Returns:
        A dictionary representation of the reference configuration, specifying
        atom positions, species mappings, and structural configuration bounds.
    """

    atoms, actual_ncells_per_side = make_fcc_template(
        ncells_per_side=ncells_per_side,
        species_list=species_list,
    )
    randomize_species(atoms, species_list, seed=seed)

    reference_config = {
        "species_list": list(species_list),
        "species_label": _species_label(species_list),
        "symbols": atoms.get_chemical_symbols(),
        "scaled_positions": atoms.get_scaled_positions().copy(),
        "ncells_per_side": int(actual_ncells_per_side),
        "seed": int(seed),
    }

    return reference_config


def generate_fcc_compute_energy_from_reference(
    model: str,
    reference_config: dict,
    alat: float,
) -> Union[Tuple[float, int], None]:
    """
    Return (total_energy, ncells_per_side) for the same FCC reference
    configuration scaled to the requested lattice constant.

    Args:
        model: Name or path of the machine learning potential model
        reference_config: The baseline crystal structure configuration used as the
          scaling template.
        alat: The target lattice parameter 'a' (in Angstroms) to scale the
          reference cell configuration to.

    Returns:
        A tuple containing the calculated total potential energy (float) and the
        number of unit cells per side (int) for the scaled system, or None if
        the energy calculation fails or encounters invalid state constraints.
    """

    symbols = reference_config["symbols"]
    scaled_positions = np.asarray(reference_config["scaled_positions"], dtype=float)
    ncells_per_side = int(reference_config["ncells_per_side"])

    atoms = Atoms(
        symbols=symbols,
        scaled_positions=scaled_positions,
        cell=[ncells_per_side * float(alat)] * 3,
        pbc=True,
    )

    calc = KIM(model)
    atoms.calc = calc

    try:
        pe = atoms.get_potential_energy()
        return float(pe), ncells_per_side

    except Exception as e:
        logger.info(
            "generate_fcc_compute_energy_from_reference exception "
            f"species={reference_config.get('species_label')} "
            f"alat={alat}:\n{e}"
        )
        return None

    finally:
        try:
            atoms.calc = None
        except Exception:
            pass

        try:
            if calc is not None and hasattr(calc, "clean"):
                calc.clean()
        except Exception:
            pass

        try:
            if calc is not None and hasattr(calc, "__del__"):
                calc.__del__()
        except Exception:
            pass

        try:
            del calc
        except Exception:
            pass

        try:
            del atoms
        except Exception:
            pass


def _coarse_scan_worker(
    model: str,
    reference_config: dict,
    a_start: float,
    a_stop: float,
    del_a: float,
    result_queue,
):
    """
    Child process for reverse coarse scan.

    Scans from a_start down to a_stop. Each successful point is sent to the
    parent immediately. If this child segfaults, the parent keeps all points
    already sent.

    Args:
        model: Name or path of the machine learning potential model
        reference_config: The baseline crystal structure configuration used as the
          starting template.
        a_start: Starting lattice parameter 'a' (in Angstroms) for the reverse
          scan.
        a_stop: Ending lattice parameter 'a' (in Angstroms) for the reverse scan.
        del_a: Step size decrement used to step through the lattice parameters.
        result_queue: A multiprocessing.Queue object used to stream successful
          energy-lattice points back to the parent process immediately upon
          evaluation.
    """

    if del_a <= 0.0:
        raise ValueError("del_a must be positive")

    if a_start < a_stop:
        raise ValueError("reverse coarse scan requires a_start >= a_stop")

    n_steps = int(math.floor((a_start - a_stop) / del_a + 1.0e-9))
    num_atoms = len(reference_config["symbols"])

    try:
        for j in range(n_steps + 1):
            alat = _round_alat(a_start - j * del_a)

            energy_config = generate_fcc_compute_energy_from_reference(
                model=model,
                reference_config=reference_config,
                alat=alat,
            )

            if energy_config is None:
                result_queue.put(
                    {
                        "type": "point",
                        "ok": False,
                        "alat": float(alat),
                        "energy_total": None,
                        "energy_per_atom": None,
                        "ncells": None,
                    }
                )
                continue

            energy_total, ncells = energy_config
            energy_per_atom = float(energy_total) / float(num_atoms)

            result_queue.put(
                {
                    "type": "point",
                    "ok": True,
                    "alat": float(alat),
                    "energy_total": float(energy_total),
                    "energy_per_atom": float(energy_per_atom),
                    "ncells": int(ncells),
                }
            )

        result_queue.put({"type": "done", "ok": True})

    except Exception:
        result_queue.put(
            {
                "type": "done",
                "ok": False,
                "error": traceback.format_exc(),
            }
        )


def _format_child_exit_reason(exitcode: Union[int, None]) -> str:
    if exitcode is None:
        return "child_exitcode_unknown"

    if exitcode == 0:
        return "completed_without_done_message"

    if exitcode < 0:
        sig = -exitcode
        try:
            sig_name = signal.Signals(sig).name
        except Exception:
            sig_name = f"signal_{sig}"
        return f"child_died_with_{sig_name}"

    return f"child_exited_with_code_{exitcode}"


def coarse_scan_reverse_safe(
    model: str,
    reference_config: dict,
    a_start: float = 12.0,
    a_stop: float = 1.5,
    del_a: float = 0.1,
    timeout: float = 600.0,
) -> dict:
    """
    Run a reverse coarse scan in one child process.

    If the child crashes or segfaults, do not restart it.
    Return all successful energy-alat points collected before the crash.

    Args:
        model: Name or path of the machine learning potential model
        reference_config: The baseline crystal structure configuration used as the
          starting template.
        a_start: Starting lattice parameter 'a' (in Angstroms) for the scan.
          Defaults to 12.0.
        a_stop: Ending lattice parameter 'a' (in Angstroms) for the scan.
          Defaults to 1.5.
        del_a: Absolute step size decrement for the lattice parameter grid points.
          Defaults to 0.1.
        timeout: Maximum execution time in seconds for the entire subprocess
          run before forcing termination. Defaults to 600.0.

    Returns:
        A dictionary containing partial or complete results up until process exit,
        including successful lattice parameters ('alats') and matched energy steps.
    """

    ctx = mp.get_context("spawn")
    result_queue = ctx.Queue()

    proc = ctx.Process(
        target=_coarse_scan_worker,
        args=(
            model,
            reference_config,
            float(a_start),
            float(a_stop),
            float(del_a),
            result_queue,
        ),
    )

    alats = []
    energies_per_atom = []
    energy_total = []
    ncells = []
    failed_alats = []

    child_done = False
    child_error = None
    child_exitcode = None
    stopped_reason = None

    def pop_from_result_queue():
        nonlocal child_done, child_error

        while True:
            try:
                result = result_queue.get_nowait()
            except queue_module.Empty:
                break

            if not isinstance(result, dict):
                continue

            if result.get("type") == "point":
                if result.get("ok", False):
                    alats.append(float(result["alat"]))
                    energies_per_atom.append(float(result["energy_per_atom"]))
                    energy_total.append(float(result["energy_total"]))
                    ncells.append(int(result["ncells"]))
                else:
                    failed_alats.append(float(result["alat"]))

            elif result.get("type") == "done":
                child_done = True
                if not result.get("ok", False):
                    child_error = result.get("error")
                break

    proc.start()
    deadline = time.monotonic() + float(timeout)

    try:
        while True:
            pop_from_result_queue()

            if child_done:
                stopped_reason = (
                    "completed" if child_error is None else "python_exception"
                )
                child_exitcode = proc.exitcode
                break

            child_exitcode = proc.exitcode

            if child_exitcode is not None:
                stopped_reason = _format_child_exit_reason(child_exitcode)

                # One last drain after child has exited.
                pop_from_result_queue()
                break

            if time.monotonic() > deadline:
                stopped_reason = "timeout"

                proc.terminate()
                proc.join(5.0)

                if proc.is_alive():
                    proc.kill()
                    proc.join()

                child_exitcode = proc.exitcode

                # Drain anything already sent before termination.
                pop_from_result_queue()
                break

            time.sleep(0.05)

    finally:
        if proc.is_alive():
            proc.terminate()
            proc.join(5.0)

            if proc.is_alive():
                proc.kill()
                proc.join()

        try:
            proc.close()
        except Exception:
            pass

        try:
            result_queue.close()
            result_queue.join_thread()
        except Exception:
            pass

    return {
        "alats": alats,
        "energies_per_atom": energies_per_atom,
        "energy_total": energy_total,
        "ncells": ncells,
        "failed_alats": failed_alats,
        "a_start": float(a_start),
        "a_stop": float(a_stop),
        "del_a": float(del_a),
        "scan_direction": "reverse",
        "stopped_reason": stopped_reason,
        "child_exitcode": child_exitcode,
        "child_error": child_error,
        "reference_config_summary": {
            "species_list": reference_config["species_list"],
            "species_label": reference_config["species_label"],
            "ncells_per_side": int(reference_config["ncells_per_side"]),
            "natoms": len(reference_config["symbols"]),
            "seed": int(reference_config["seed"]),
        },
    }


def _local_minima_indices(energies_per_atom: list[float]) -> list[int]:
    """
    Find the array indices of local energy minima using peak detection.

    Inverts the energy values to transform local minima into local maxima peaks,
    then locates them using SciPy's peak finder while discarding any non-finite
    results.

    Args:
        energies_per_atom: A sequence of calculated energy values per atom over
          the scan grid.

    Returns:
        A list of integer indices pinpointing the locations of valid local minima
        within the input sequence.
    """

    y = np.asarray(energies_per_atom, dtype=float)

    if len(y) == 0 or np.sum(np.isfinite(y)) == 0:
        return []

    indices, _ = find_peaks(-y)

    return [int(i) for i in indices if np.isfinite(y[i])]


def _energy_is_within_bounds(
    energy_per_atom: float,
    energy_bound: tuple[float, float],
) -> bool:
    """
    Check if the magnitude of the given energy value falls within specific bounds.

    Ensures that the parsed energy value is a finite numerical value before
    comparing its absolute magnitude against the minimum and maximum limits.

    Args:
        energy_per_atom: The calculated energy value per atom to evaluate.
        energy_bound: A tuple containing the (minimum, maximum) acceptable absolute
          energy thresholds.

    Returns:
        True if the absolute energy magnitude is finite and within the bounds,
        otherwise False.
    """

    energy_per_atom = float(energy_per_atom)

    if not np.isfinite(energy_per_atom):
        return False

    energy_magnitude = abs(energy_per_atom)

    return float(energy_bound[0]) <= energy_magnitude <= float(energy_bound[1])


def _starting_points_from_scan(
    scan: dict,
    max_starting_points: int,
    energy_bound: tuple[float, float],
) -> list[dict]:
    """
    Identify and sort viable starting points for optimization from a coarse scan.

    Finds local energy minima from scan data, filters them by energy limits,
    and ranks them from lowest to highest energy to prioritize the most stable
    configurations.

    Args:
        scan: Dictionary containing lists of lattice parameters ('alats') and
          corresponding 'energies_per_atom' from the coarse sweep.
        max_starting_points: Maximum number of candidate minima to return after
          sorting.
        energy_bound: Valid range (min, max) for acceptable configuration
          energies per atom.

    Returns:
        A list of dictionaries, where each item contains the original index,
        lattice parameter ('alat'), and energy value of a valid starting point.
    """

    minima_indices = _local_minima_indices(scan["energies_per_atom"])

    minima_indices = [
        i
        for i in minima_indices
        if _energy_is_within_bounds(
            scan["energies_per_atom"][i],
            energy_bound,
        )
    ]

    minima_indices = sorted(
        minima_indices,
        key=lambda i: scan["energies_per_atom"][i],
    )

    minima_indices = minima_indices[:max_starting_points]

    return [
        {
            "index": int(i),
            "alat": float(scan["alats"][i]),
            "energy_per_atom": float(scan["energies_per_atom"][i]),
        }
        for i in minima_indices
    ]


def _nelder_mead_worker(
    model_name: str,
    reference_config: dict,
    start_alat: float,
    a_min: float,
    a_max: float,
    energy_bound: tuple[float, float],
    queue,
):
    """
    Run Nelder-Mead in a child process using the same reference FCC
    configuration as the coarse scan.

    Args:
        model_name: Name or path of the machine learning potential model
        reference_config: The baseline crystal structure configuration used as the
          starting template.
        start_alat: Initial lattice parameter 'a' value (in Angstroms) to begin the
          optimization.
        a_min: Minimum allowable bounding value for the lattice parameter 'a'.
        a_max: Maximum allowable bounding value for the lattice parameter 'a'.
        energy_bound: Valid range (min, max) for acceptable configuration
          energies.
        queue: A multiprocessing.Queue object used to safely pass the optimization
          results back to the parent process.
    """

    try:
        evaluations = []
        num_atoms = len(reference_config["symbols"])

        def objective(x):
            alat = float(np.ravel(x)[0])

            if not np.isfinite(alat) or alat < a_min or alat > a_max:
                return 1.0e100

            val = generate_fcc_compute_energy_from_reference(
                model=model_name,
                reference_config=reference_config,
                alat=alat,
            )

            if val is None:
                return 1.0e100

            energy_total, ncells = val
            energy_per_atom = float(energy_total) / float(num_atoms)

            if not _energy_is_within_bounds(energy_per_atom, energy_bound):
                return 1.0e100

            evaluations.append(
                {
                    "alat": float(alat),
                    "energy_total": float(energy_total),
                    "energy_per_atom": float(energy_per_atom),
                    "ncells": int(ncells),
                }
            )

            return float(energy_per_atom)

        result = minimize(
            objective,
            x0=np.array([float(start_alat)]),
            method="Nelder-Mead",
            options={
                "xatol": 1.0e-4,
                "fatol": 1.0e-8,
                "maxiter": 80,
                "maxfev": 160,
                "disp": False,
            },
        )

        if len(evaluations) == 0:
            queue.put(
                {
                    "ok": False,
                    "status": "no_valid_evaluations",
                    "message": str(result.message),
                    "start_alat": float(start_alat),
                    "evaluations": evaluations,
                }
            )
            return

        best_eval = min(evaluations, key=lambda row: row["energy_per_atom"])

        if not result.success:
            queue.put(
                {
                    "ok": False,
                    "status": "optimizer_unsuccessful",
                    "message": str(result.message),
                    "start_alat": float(start_alat),
                    "best_alat": float(best_eval["alat"]),
                    "best_energy_per_atom": float(best_eval["energy_per_atom"]),
                    "evaluations": evaluations,
                }
            )
            return

        queue.put(
            {
                "ok": True,
                "status": "success",
                "message": str(result.message),
                "start_alat": float(start_alat),
                "good_alat": float(best_eval["alat"]),
                "good_energy_total": float(best_eval["energy_total"]),
                "good_energy_per_atom": float(best_eval["energy_per_atom"]),
                "good_ncells": int(best_eval["ncells"]),
                "optimizer_x": float(np.ravel(result.x)[0]),
                "optimizer_fun": float(result.fun),
                "nfev": int(result.nfev),
                "nit": int(result.nit),
                "evaluations": evaluations,
            }
        )

    except Exception:
        queue.put(
            {
                "ok": False,
                "status": "python_exception",
                "message": traceback.format_exc(),
                "start_alat": float(start_alat),
                "evaluations": [],
            }
        )


def scipy_nelder_mead_safe(
    model: str,
    reference_config: dict,
    start_alat: float,
    a_min: float,
    a_max: float,
    energy_bound: tuple[float, float],
    timeout: float = 600.0,
) -> dict:
    """
    Run Nelder-Mead in a child process.

    This isolates crashes/segfaults from the parent process.

    Args:
        model: Name or path of the machine learning potential model
        reference_config: The baseline crystal structure configuration used as the
          starting template.
        start_alat: Initial lattice parameter 'a' value (in Angstroms) to begin the
          optimization.
        a_min: Minimum allowable bounding value for the lattice parameter 'a'.
        a_max: Maximum allowable bounding value for the lattice parameter 'a'.
        energy_bound: Valid range (min, max) for acceptable configuration
          energies.
        timeout: Maximum execution time in seconds before terminating the child
          process worker. Defaults to 600.0.

    Returns:
        A dictionary containing the optimization outcome status, evaluations, and
        the refined equilibrium lattice coordinates if successful.
    """

    def fail_result(status: str, message: str) -> dict:
        return {
            "ok": False,
            "status": status,
            "message": message,
            "start_alat": float(start_alat),
            "evaluations": [],
        }

    ctx = mp.get_context("spawn")
    result_queue = ctx.Queue(maxsize=1)

    proc = ctx.Process(
        target=_nelder_mead_worker,
        args=(
            model,
            reference_config,
            start_alat,
            a_min,
            a_max,
            energy_bound,
            result_queue,
        ),
    )

    try:
        proc.start()
        proc.join(timeout)

        if proc.is_alive():
            proc.terminate()
            proc.join(5.0)

            if proc.is_alive():
                proc.kill()
                proc.join()

            logger.info(
                "scipy_nelder_mead_safe TIMEOUT for "
                f"{model} {reference_config['species_list']} "
                f"start_alat={start_alat}"
            )

            return fail_result(
                "timeout",
                f"Nelder-Mead timed out after {timeout} seconds",
            )

        if proc.exitcode != 0:
            message = _format_child_exit_reason(proc.exitcode)

            logger.info(
                "scipy_nelder_mead_safe CRASH/EXIT for "
                f"{model} {reference_config['species_list']} "
                f"start_alat={start_alat}: {message}"
            )

            return fail_result("crash", message)

        try:
            result = result_queue.get(timeout=1.0)
        except queue_module.Empty:
            return fail_result(
                "no_result",
                "child exited but returned no result",
            )

        if not isinstance(result, dict):
            return fail_result(
                "invalid_result",
                f"child returned invalid result type {type(result).__name__}",
            )

        if "ok" not in result:
            return fail_result(
                "invalid_result",
                "child result dictionary is missing required key 'ok'",
            )

        if not result.get("ok", False):
            result.setdefault("status", "worker_error")
            result.setdefault("message", "Nelder-Mead worker returned ok=False")
            result.setdefault("start_alat", float(start_alat))
            result.setdefault("evaluations", [])
            return result

        result.setdefault("status", "success")
        result.setdefault("message", "")
        result.setdefault("start_alat", float(start_alat))
        result.setdefault("evaluations", [])

        return result

    finally:
        if proc.is_alive():
            proc.terminate()
            proc.join(5.0)

            if proc.is_alive():
                proc.kill()
                proc.join()

        try:
            result_queue.close()
            result_queue.join_thread()
        except Exception:
            pass

        try:
            proc.close()
        except Exception:
            pass


def _failure_config_result(
    species_list: list[str],
    configuration_type: str,
    reference_config: Union[dict, None],
    scan: Union[dict, None],
    starting_points: list[dict],
    attempts: list[dict],
    reason: str,
    bounds: dict,
) -> dict:
    """Generate a uniform dictionary summary for a failed configuration search.

    Populates default error indicators like negative dimensions or None values for
    missing measurements while preserving logs, boundaries, and historical attempts
    for later diagnosis.

    Args:
        species_list: Chemical symbols of the elements included in the system
        configuration_type: A label indicating the structural type
        reference_config: The generated baseline crystal structure, or None if
          generation failed.
        scan: Historical trace of the initial energy-versus-lattice coarse scan,
          or None if skipped.
        starting_points: Local energy minima detected during the coarse scan to
          be used for optimization.
        attempts: List of optimization traces produced by downstream Nelder-Mead
          subprocess workers.
        reason: A string identifier explaining why the search workflow halted.
        bounds: Dictionary defining coordinate search window ranges and step delta
          sizes.

    Returns:
        A dictionary containing error statuses, tracking identifiers, and any
        partial execution traces collected before the failure occurred.
    """

    result = {
        "ok": False,
        "status": "failed",
        "failure_reason": reason,
        "configuration_type": configuration_type,
        "species_list": list(species_list),
        "species_label": _species_label(species_list),
        "good_alat": -1.0,
        "good_energy_per_atom": None,
        "good_energy_total": None,
        "good_ncells": -1,
        "bounds": bounds,
        "coarse_scan": scan,
        "coarse_minima": starting_points,
        "nelder_mead_attempts": attempts,
        "search_strategy": "reverse_coarse_scan_plus_first_successful_nelder_mead",
        "all_alats": [],
        "all_energies_per_atom": [],
    }

    if reference_config is not None:
        result["reference_config"] = reference_config

    return result


def _equilibrate_one_config(
    model: str,
    species_list: list[str],
    configuration_type: str,
    energy_bound: tuple[float, float],
    coarse_del_a: float,
    coarse_timeout: float,
    nelder_mead_timeout: float,
    max_starting_points: int,
    ncells_per_side: int = DEFAULT_FCC_NCELLS_PER_SIDE,
    seed: int = 13,
    coarse_a_start: float = 12.0,
    coarse_a_stop: float = 1.5,
) -> dict:
    """Equilibrate a single FCC atomic configuration using a sequential protocol.

    Executes a high-level orchestration that generates a reference crystal,
    scans lattice constants to find energy minima, and applies a Nelder-Mead
    simplex search to optimize structure bounds.

    Args:
        model: Name or path of the machine learning potential model
        species_list: Chemical symbols of the elements to include in the system
        configuration_type: A label indicating the structural type
        energy_bound: Valid range (min, max) for acceptable configuration
          energies.
        coarse_del_a: Step size for the lattice parameter 'a' during coarse
          scans.
        coarse_timeout: Maximum execution time in seconds for the coarse scan
          subprocess.
        nelder_mead_timeout: Maximum execution time in seconds for each
          Nelder-Mead optimization subprocess.
        max_starting_points: Maximum number of local minima to evaluate using
          Nelder-Mead relaxation.
        ncells_per_side: Number of unit cells along each axis to define the
          supercell size.
        seed: Random seed for reproducible structure generation and atom
          placement.
        coarse_a_start: Starting lattice parameter 'a' (in Angstroms) for the
          reverse scan.
        coarse_a_stop: Ending lattice parameter 'a' (in Angstroms) for the
          reverse scan.

    Returns:
        A dictionary containing the metadata, execution history, and details of
        either a successful equilibrium search or a recorded failure.
    """

    logger.info(
        f"_equilibrate_one_config {_species_label(species_list)} "
        f"config_type={configuration_type} "
        f"coarse_a_start {coarse_a_start}, coarse_a_stop {coarse_a_stop} "
        f"ncells_per_side {ncells_per_side} "
    )

    reference_config = make_fcc_reference_config(
        species_list=species_list,
        ncells_per_side=ncells_per_side,
        seed=seed,
    )

    scan = coarse_scan_reverse_safe(
        model=model,
        reference_config=reference_config,
        a_start=coarse_a_start,
        a_stop=coarse_a_stop,
        del_a=coarse_del_a,
        timeout=coarse_timeout,
    )

    bounds = {
        "coarse_a_start": float(coarse_a_start),
        "coarse_a_stop": float(coarse_a_stop),
        "coarse_del_a": float(coarse_del_a),
        "nelder_mead_a_min": float(min(coarse_a_start, coarse_a_stop)),
        "nelder_mead_a_max": float(max(coarse_a_start, coarse_a_stop)),
    }

    if len(scan["alats"]) == 0:
        return _failure_config_result(
            species_list,
            configuration_type,
            reference_config,
            scan,
            [],
            [],
            "coarse_scan_has_no_successful_points",
            bounds,
        )

    starting_points = _starting_points_from_scan(
        scan=scan,
        max_starting_points=max_starting_points,
        energy_bound=energy_bound,
    )

    if len(starting_points) == 0:
        reason = "no_coarse_minima_found"

        if scan.get("stopped_reason") not in (
            "completed",
            "completed_without_done_message",
        ):
            reason = f"no_coarse_minima_found_before_{scan.get('stopped_reason')}"

        return _failure_config_result(
            species_list,
            configuration_type,
            reference_config,
            scan,
            [],
            [],
            reason,
            bounds,
        )

    attempts = []

    for start_index, start in enumerate(starting_points):
        attempt = scipy_nelder_mead_safe(
            model=model,
            reference_config=reference_config,
            start_alat=start["alat"],
            a_min=bounds["nelder_mead_a_min"],
            a_max=bounds["nelder_mead_a_max"],
            energy_bound=energy_bound,
            timeout=nelder_mead_timeout,
        )

        attempts.append(attempt)

        if attempt.get("ok", False):
            evaluations = attempt.get("evaluations", [])

            return {
                "ok": True,
                "status": "success",
                "configuration_type": configuration_type,
                "species_list": list(species_list),
                "species_label": _species_label(species_list),
                "good_alat": float(attempt["good_alat"]),
                "good_energy_per_atom": float(attempt["good_energy_per_atom"]),
                "good_energy_total": float(attempt["good_energy_total"]),
                "good_ncells": int(attempt["good_ncells"]),
                "bounds": bounds,
                "reference_config": reference_config,
                "coarse_scan": scan,
                "coarse_minima": starting_points,
                "nelder_mead_attempts": attempts,
                "selected_attempt_index": int(start_index),
                "selected_policy": "first_successful_nelder_mead",
                "search_strategy": (
                    "reverse_coarse_scan_plus_first_successful_nelder_mead"
                ),
                "all_alats": [float(row["alat"]) for row in evaluations],
                "all_energies_per_atom": [
                    float(row["energy_per_atom"]) for row in evaluations
                ],
            }

    return _failure_config_result(
        species_list,
        configuration_type,
        reference_config,
        scan,
        starting_points,
        attempts,
        "all_nelder_mead_attempts_failed",
        bounds,
    )


def _finite_positive_mean(values: list[float]) -> float:
    """Calculate the mean of finite, positive numbers in a list.

    Filters out non-finite numbers (NaN, Inf) and values less than or equal to
    zero before computing the average.

    Args:
        values: A list of floats to be filtered and averaged.

    Returns:
        The arithmetic mean of the remaining valid values as a float, or -1.0
        if no values survive the filtering process.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    arr = arr[arr > 0.0]

    if len(arr) == 0:
        return -1.0

    return float(np.mean(arr))


def _minimal_success_result(full_result: dict) -> dict:
    """Extract a minimal summary from a successful full result dictionary.

    Args:
        full_result: The complete dictionary of results containing metadata,
          scan details, and final optimization outputs.

    Returns:
        A dictionary containing parsed information about the equilibrium
        configuration, including energy values, cell parameters, and structural
        metadata.
    """

    final_result = full_result["final_result"]

    return {
        "ok": True,
        "status": "success",
        "model": full_result["model"],
        "species_list": full_result["species_list"],
        "species_label": full_result["species_label"],
        "configuration_type": final_result["configuration_type"],
        "equilibrium_alat": float(final_result["good_alat"]),
        "approx_mixed_equilibrium_alat": float(
            full_result["approx_mixed_equilibrium_alat"]
        ),
        "mono_species_equilibrium_alats": full_result["mono_species_equilibrium_alats"],
        "ncells_per_side": int(final_result["good_ncells"]),
        "natoms": fcc_atoms_in_supercell(int(final_result["good_ncells"])),
        "energy_per_atom": float(final_result["good_energy_per_atom"]),
        "energy_total": float(final_result["good_energy_total"]),
        "method": full_result["method"],
    }


def find_equilibrium_config_FCC(
    model: str,
    species_list: list[str],
    energy_bound: tuple[float, float] = (5.0e-2, 5.0e2),
    coarse_del_a: float = 0.1,
    mixed_coarse_del_a: float = 0.1,
    coarse_timeout: float = 600.0,
    nelder_mead_timeout: float = 600.0,
    max_starting_points: int = 6,
    ncells_per_side: int = DEFAULT_FCC_NCELLS_PER_SIDE,
    seed: int = 13,
    coarse_a_start: float = 12.0,
    coarse_a_stop: float = 1.5,
) -> dict:
    """
    Find an FCC equilibrium configuration.

    Protocol:
      1. For each mono species, create one fixed FCC reference configuration.
      2. Run one reverse coarse scan subprocess from coarse_a_start to
         coarse_a_stop.
      3. Keep all energy-alat points collected before any child crash.
      4. Find local minima using scipy.signal.find_peaks(-E).
      5. Run Nelder-Mead subprocesses from minima sorted by coarse energy.
      6. Return after the first successful Nelder-Mead result.
      7. For mixed species, repeat the same protocol using one fixed mixed
         reference configuration.

    Args:
        model: Name or path of the machine learning potential model.
        species_list: Chemical symbols of the elements to include in the system
        energy_bound: Valid range (min, max) for configuration energies.
        coarse_del_a: Step size for the lattice parameter 'a' for monospecies
          coarse scans.
        mixed_coarse_del_a: Step size for the lattice parameter 'a' for
          mixed-species coarse scans.
        coarse_timeout: Maximum execution time in seconds for the coarse scan
          subprocess.
        nelder_mead_timeout: Maximum execution time in seconds for each
          Nelder-Mead optimization subprocess.
        max_starting_points: Maximum number of local minima to evaluate using
          Nelder-Mead relaxation.
        ncells_per_side: Number of unit cells along each axis to define the
          supercell size.
        seed: Random seed for reproducible structure generation and atom
          placement.
        coarse_a_start: Starting lattice parameter 'a' (in Angstroms) for the
          reverse scan.
        coarse_a_stop: Ending lattice parameter 'a' (in Angstroms) for the
          reverse scan.

    Returns:
        A dictionary containing the successful equilibrium configuration results
        and metadata.

    """

    mono_results = []

    for species in species_list:
        mono = _equilibrate_one_config(
            model=model,
            species_list=[species],
            configuration_type="mono_species_fcc",
            energy_bound=energy_bound,
            coarse_del_a=coarse_del_a,
            coarse_timeout=coarse_timeout,
            nelder_mead_timeout=nelder_mead_timeout,
            max_starting_points=max_starting_points,
            ncells_per_side=1,
            seed=seed,
            coarse_a_start=coarse_a_start,
            coarse_a_stop=coarse_a_stop,
        )

        mono_results.append(mono)

    mono_good_alats = [row.get("good_alat", -1.0) for row in mono_results]
    approx_mixed_equilibrium_alat = _finite_positive_mean(mono_good_alats)

    if len(species_list) == 1:
        mixed_result = None
        final_result = mono_results[0]
        results = mono_results

    else:
        mixed_result = _equilibrate_one_config(
            model=model,
            species_list=species_list,
            configuration_type="mixed_species_fcc",
            energy_bound=energy_bound,
            coarse_del_a=mixed_coarse_del_a,
            coarse_timeout=coarse_timeout,
            nelder_mead_timeout=nelder_mead_timeout,
            max_starting_points=max_starting_points,
            ncells_per_side=1,
            seed=seed,
            coarse_a_start=coarse_a_start,
            coarse_a_stop=coarse_a_stop,
        )

        mixed_result["approx_mixed_equilibrium_alat"] = approx_mixed_equilibrium_alat

        final_result = mixed_result
        results = mono_results + [mixed_result]

    full_result = {
        "ok": bool(final_result and final_result.get("ok", False)),
        "status": (
            "success" if final_result and final_result.get("ok", False) else "failed"
        ),
        "model": model,
        "species_list": list(species_list),
        "species_label": _species_label(species_list),
        "ncells_per_side": ncells_per_side,
        "method": "reverse_coarse_scan_plus_first_successful_nelder_mead",
        "mono_species_equilibrium_alats": {
            row["species_label"]: row.get("good_alat", -1.0) for row in mono_results
        },
        "approx_mixed_equilibrium_alat": approx_mixed_equilibrium_alat,
        "mono_species_results": mono_results,
        "mixed_species_result": mixed_result,
        "final_result": final_result,
        "equilibrium_alat": (
            final_result.get("good_alat", -1.0) if final_result else -1.0
        ),
        "results": results,
    }

    if full_result["ok"]:
        return _minimal_success_result(full_result)

    return full_result


################################################################################

# If called directly, do nothing
if __name__ == "__main__":
    pass
