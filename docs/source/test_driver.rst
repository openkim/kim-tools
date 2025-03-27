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


You must create a class named ``TestDriver`` inheriting from :class:`~kim_tools.SingleCrystalTestDriver`.
In your ``TestDriver`` class, you must overload the function :func:`~kim_tools.KIMTestDriver._calculate`. 
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
Test Driver in order for your Driver to work in the OpenKIM pipeline: :doc:`tutorial_pipeline`. Before this, you will likely wish
to debug your code by invoking it directly in Python. See how to do that here: :doc:`tutorial_debug`

**Please read the code and the comments carefully, as they explain many important aspects of writing your own Test Driver.** You should understand
the usage of the following functions. Click the links below for more information on each:

- ``self.``:func:`~kim_tools.SingleCrystalTestDriver._get_atoms`:
  returns an :class:`~ase.Atoms` object representing a primitive unit cell of the crystal as a starting point for your calculations.

- ``self.``:func:`~kim_tools.SingleCrystalTestDriver._update_nominal_parameter_values`: 
  if your Test Driver changes the crystal structure, pass this a primitive unit cell of the crystal to update the nominal crystal structure.
  
  .. note::

    If the symmetry of your crystal has changed in the process of your simulation, this function will raise an error.
    This is intended and it is expected that you do not handle this exception. It is our convention that if the 
    crystal undergoes a symmetry-changing phase transition, the result is invalid and the Test Driver should exit with an error.

- ``self.``:func:`~kim_tools.SingleCrystalTestDriver._verify_unchanged_symmetry`:
  You may also find it useful to check for a symmetry change without changing the stored description of the crystal. This is useful, for example,
  if your Test Driver is changing temperature or stress around some reference value, and you wish to make sure that the state changes have not
  induced a phase transition. This function takes an :class:`~ase.Atoms` object and returns ``True`` if the symmetry has not changed, otherwise ``False``.

- ``self.``:func:`~kim_tools.SingleCrystalTestDriver._add_property_instance_and_common_crystal_genome_keys`:
  Use this to initialize an Instance of the KIM Property(s) you defined in :doc:`tutorial_property`. It will automatically be populated with the keys
  describing the nominal state of the crystal.

- ``self.``:func:`~kim_tools.SingleCrystalTestDriver._add_key_to_current_property_instance`:
  Use this to add additional keys to the last Property Instance you created.

- ``self.``:func:`~kim_tools.SingleCrystalTestDriver._get_temperature`

- ``self.``:func:`~kim_tools.SingleCrystalTestDriver._get_cell_cauchy_stress`

- ``self.``:func:`~kim_tools.SingleCrystalTestDriver._get_nominal_crystal_structure_npt`:
  If you need fine-grained access to the symmetry-reduced AFLOW prototype designation of the crystal


.. _doc.example_test_driver:

Example ``test_driver.py``
==========================

.. literalinclude:: ../../examples/CrystalGenomeASEExample__TD_000000654321_000/test_driver/test_driver.py
    :language: Python

Functionality not covered in the above example
==============================================

- ``self.``:func:`~kim_tools.KIMTestDriver._calc`: 
  This gives you access to the ASE calculator object, if you are building a separate :class:`~ase.Atoms` object and need to attach a calculator to it.

- ``self.``:attr:`~kim_tools.KIMTestDriver.kim_model_name`:
  The KIM Model Name (if present). You should only use this if exporting data to a non-ASE simulator (see below).

- ``self.``:func:`~kim_tools.KIMTestDriver._add_file_to_current_property_instance`:
  This adds a "file" type key to the current Property Instance, for example the "restart-file" key mentioned in the note below. It will automatically
  be numbered according to the current Property Instance and moved to the ``output`` directory if it is not already there (e.g. ``restart.dump``
  will be automatically moved to the path ``output/restart-1.dump`` and reported in the property accordingly.)

.. note::

  If you are running an MD simulation, the structure you report should be time-averaged and likely averaged over the supercell folded back into the unit cell. 
  This will give you more robust averages. Additionally, if you are performing an NPT simulation, you may as well write an instance of the 
  `crystal-structure-npt <https://openkim.org/properties/show/crystal-structure-npt>`_ property for future re-use. Note the optional ``restart-file`` key. 
  It is recommended that you save a restart file (for example, ``restart.dump``). You can then add it using ``self.``:func:`~kim_tools.KIMTestDriver._add_file_to_current_property_instance`.

  .. code-block:: Python

    self._add_property_instance_and_common_crystal_genome_keys("crystal-structure-npt",write_temp=True,write_stress=True)
    self._add_file_to_current_property_instance("restart-file","restart.dump")  

.. todo::
    
  This should be an integrated part of ``kim-tools``

LAMMPS and other non-ASE simulators
-----------------------------------
.. note::

  Non-ASE calculations require extra steps. 
  For example, one way to run a LAMMPS simulation in this framework is to export your atomic configuration using :func:`ase.io.write`, create LAMMPS input file(s)
  with `kim commands <https://docs.lammps.org/kim_commands.html>`_ using the KIM model stored in the base class' attribute ``self.``:attr:`~kim_tools.KIMTestDriver.kim_model_name`,
  run your simulation(s), and read the configuration back in using :func:`ase.io.read` (for example, to re-detect the changed crystal structure).

.. todo::
    
  This should be an integrated part of ``kim-tools``
