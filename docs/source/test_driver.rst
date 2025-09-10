=======================
Main Computational Code
=======================

.. contents:: Table of Contents

The next step is to write the main computational code. It must be contained in a Python file named ``test_driver/test_driver.py``,
although you may include as many other Python files as you wish for utility functions.

.. note::

  Your ``test_driver`` directory should be structured as a Python package. This means that it contains a (probably empty) ``__init__.py``,
  and you should use relative imports to import any other Python files from your ``test_driver.py``. For example, if you want to put
  some function ``my_helper_function`` in a Python file ``helper_functions.py``, you shoud structure your ``test_driver`` directory like this:

  ::

    test_driver/
    ├─ __init__.py
    ├─ test_driver.py
    ├─ helper_functions.py

  Then, in your ``test_driver.py`` you would import ``my_helper_function`` like this:

  .. code-block:: Python

    from .helper_functions import my_helper_function


You must create a class named ``TestDriver`` inheriting from :class:`~kim_tools.test_driver.core.SingleCrystalTestDriver`.
In your ``TestDriver`` class, you must overload the function :func:`~kim_tools.test_driver.core.KIMTestDriver._calculate`.
Besides ``self``, the function must also accept ``**kwargs``. Before ``**kwargs``, you may add any additional arguments that
you wish users or the OpenKIM Pipeline to be able to vary.

.. note::

  Temperature and stress are commonly used, so they do not need to be added as additional arguments. If needed,
  your Test Driver will automatically accept arguments ``temperature_K`` and ``cell_cauchy_stress_eV_angstrom3``. See the
  :ref:`doc.example_test_driver` for an example of how to access them from inside your ``_calculate`` function.

.. todo::

  Replace ``temperature_K`` and ``cell_cauchy_stress_eV_angstrom3`` with ``temperature`` and ``cell_cauchy_stress``, and
  allow the user to specify units?

An example ``test_driver.py`` from |example_url| is shown below. The ``_calculate`` function computes the energy vs. volume
curve for isotropic expansion and compression of a crystal at zero temperature. You can use it as a starting point for your Test Driver.
Additional documentation about functionality not demonstrated in the example can be found at the bottom of this page.

Once ``test_driver.py`` is written, you will have to make some small modifications to the auxiliary files provided in the example
Test Driver in order for your Driver to work in the OpenKIM pipeline and the kimvv package: :doc:`tutorial_pipeline`, :doc:`tutorial_package`. Before this, you will likely wish
to debug your code by invoking it directly in Python. See how to do that here: :doc:`tutorial_debug`

**Please read the code and the comments carefully, as they explain many important aspects of writing your own Test Driver.** You should understand
the usage of the following functions. Click the links below for more information on each:

- ``self.``:func:`~kim_tools.test_driver.core.SingleCrystalTestDriver._get_atoms`:
  returns an :class:`~ase.Atoms` object representing a primitive unit cell of the crystal as a starting point for your calculations.

- ``self.``:func:`~kim_tools.test_driver.core.SingleCrystalTestDriver._update_nominal_parameter_values`:
  if your Test Driver changes the crystal structure, pass this a primitive unit cell of the crystal to update the nominal crystal structure.

  .. note::

    If the symmetry of your crystal has changed in the process of your simulation, this function will raise an error.
    This is intended and it is expected that you do not handle this exception. It is our convention that if the
    crystal undergoes a symmetry-changing phase transition, the result is invalid and the Test Driver should exit with an error.

- ``self.``:func:`~kim_tools.test_driver.core.SingleCrystalTestDriver._verify_unchanged_symmetry`:
  You may also find it useful to check for a symmetry change without changing the stored description of the crystal. This is useful, for example,
  if your Test Driver is changing temperature or stress around some reference value, and you wish to make sure that the state changes have not
  induced a phase transition. This function takes an :class:`~ase.Atoms` object and returns ``True`` if the symmetry has not changed, otherwise ``False``.

- ``self.``:func:`~kim_tools.test_driver.core.SingleCrystalTestDriver._add_property_instance_and_common_crystal_genome_keys`:
  Use this to initialize an Instance of the KIM Property(s) you defined in :doc:`tutorial_property`. It will automatically be populated with the keys
  describing the nominal state of the crystal.

- ``self.``:func:`~kim_tools.test_driver.core.KIMTestDriver._add_key_to_current_property_instance`:
  Use this to add additional keys to the last Property Instance you created.

- ``self.``:func:`~kim_tools.test_driver.core.SingleCrystalTestDriver._get_temperature`

- ``self.``:func:`~kim_tools.test_driver.core.SingleCrystalTestDriver._get_cell_cauchy_stress`

- ``self.``:func:`~kim_tools.test_driver.core.SingleCrystalTestDriver._get_nominal_crystal_structure_npt`:
  If you need fine-grained access to the symmetry-reduced AFLOW prototype designation of the crystal

- ``self.``:func:`~kim_tools.test_driver.core.KIMTestDriver._add_file_to_current_property_instance`:
  For adding keys with the "file" type to your Property Instance, e.g. restart files or images


.. _doc.example_test_driver:

Example ``test_driver.py``
==========================

.. literalinclude:: ../../examples/CrystalGenomeASEExample__TD_000000654321_000/test_driver/test_driver.py
    :language: Python

Functionality not covered in the above example
==============================================

- ``self.``:func:`~kim_tools.test_driver.core.KIMTestDriver._calc`:
  This gives you access to the ASE calculator object, if you are building a separate :class:`~ase.Atoms` object and need to attach a calculator to it.

- ``self.``:attr:`~kim_tools.test_driver.core.KIMTestDriver.kim_model_name`:
  The KIM Model Name (if present). You should only use this if exporting data to a non-ASE simulator (see below).

Molecular Dynamics
------------------

If you are running an MD simulation, the structure you report should be time-averaged and likely averaged over the supercell folded back into the unit cell.
This will give you more robust averages assuming the translational symmetry of the unit cell was not broken. At this time, the functions
:func:`kim_tools.symmetry_util.core.reduce_and_avg` and :func:`kim_tools.symmetry_util.core.kstest_reduced_distances` allow you to perform these operations,
assuming your supercell is built from contiguous repeats of the unit cell (i.e. atoms 0 to *N*-1 in the supercell are the original unit cell, atoms
*N* to 2 *N*-1 are a the original unit cell shifted by an integer multiple of the lattice vectors, and so on). See
https://github.com/openkim/kim-tools/blob/main/tests/test_symmetry_util.py for an example of how to use these functions.

Additionally, if you are performing an NPT simulation, you may as well write an instance of the
`crystal-structure-npt <https://openkim.org/properties/show/crystal-structure-npt>`_ property for future re-use. Note the optional ``restart-file`` key.
It is recommended that you save a restart file (for example, ``restart.dump``).
You can then add it using ``self.``:func:`~kim_tools.test_driver.core.KIMTestDriver._add_file_to_current_property_instance`.

.. code-block:: Python

  self._add_property_instance_and_common_crystal_genome_keys("crystal-structure-npt",write_temp=True,write_stress=True)
  self._add_file_to_current_property_instance("restart-file","restart.dump")

.. todo::

  Update :func:`kim_tools.symmetry_util.core.change_of_basis_atoms` to incorporate the functionalities of :func:`kim_tools.symmetry_util.core.reduce_and_avg`
  and :func:`kim_tools.symmetry_util.core.kstest_reduced_distances`, without the restriction on how the supercell is built.

Writting Test Drivers using LAMMPS
----------------------------------

In general, using ASE to perform the computations is preferable, but you may need access to LAMMPS functionality, for example to run a large parallel MD or static
calculation with domain decomposition. You may use LAMMPS in whichever way is most convenient for you as long as it is wrapped within the ``_calculate`` Python
method.
For example, one way to run a LAMMPS simulation in this framework is to export your atomic configuration using :func:`ase.io.write`, create LAMMPS input file(s)
with `kim commands <https://docs.lammps.org/kim_commands.html>`_ using the KIM model stored in the base class' attribute ``self.``:attr:`~kim_tools.test_driver.core.KIMTestDriver.kim_model_name`,
run your simulation(s), and read the configuration back in using :func:`ase.io.read` (for example, to re-detect the changed crystal structure).

There are some special considerations when using LAMMPS to write KIM Test Drivers. See https://github.com/openkim-hackathons/CrystalGenomeLAMMPSExample__TD_000000654322_000
for a trivial Crystal Genome LAMMPS example.

* Not all models support all LAMMPS systems of units. To work around this, use the ``unit_conversion_mode`` option for the ``kim init`` command and use the
  LAMMPS variables created by this command to convert to a single set of units. See https://docs.lammps.org/kim_commands.html#openkim-im-initialization-kim-init
  for more info.
* The format of a LAMMPS data file must be commensurate to the LAMMPS ``atom_style`` declared. Most KIM models use ``atom_style atomic``, however some use ``atom_style charge``.
  To determine the ``atom_style``, use ``self.``:func:`~kim_tools.test_driver.core.KIMTestDriver._get_supported_lammps_atom_style`. You can then pass the value returned by this
  function to :func:`ase.io.write` to generate a LAMMPS data file with the appropriate ``atom_style``.
* By default, LAMMPS will un-skew the box if it gets too tilted. This can cause :func:`~kim_tools.test_driver.core.SingleCrystalTestDriver._update_nominal_parameter_values`
  to fail. To suppress this LAMMPS behavior, use the ``flip no`` option with ``fix npt`` or ``fix deform``. See
  https://docs.lammps.org/Howto_triclinic.html#periodicity-and-tilt-factors-for-triclinic-simulation-boxes for more info.

.. todo::

  LAMMPS capabilities should be better integrated into ``kim-tools`` and better examples should be provided.
