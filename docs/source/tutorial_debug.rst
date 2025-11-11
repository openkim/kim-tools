=====================================================
Testing and Debugging Your Crystal Genome Test Driver
=====================================================

When your Test Driver is integrated into the OpenKIM Processing Pipeline, the auxiliary files supplied with your
Test Driver will be used to pass inputs and parse outputs from the calculation you have implemented in ``test_driver.py``.
This functionality is also available for testing in the KIM Developer Platform. More information about this is avalable
here: :doc:`../../tutorial_pipeline`, however it is likely that you will first wish to test your Driver by invoking it directly.

Examples for how to do that are found in the script ``run.py`` from |example_url|. Use this as a starting point for running your
own Test Driver, changing the passed arguments and the properties printed out. To see a rendered version of that file explaining its functionality,
follow the link in the :ref:`doc.example_script` section below.

A practical guide for testing your Driver on a variety of crystal structures is found at the bottom of the page in the section ':ref:`doc.curated_tests`'.
Everything you need to run the example script is containerized in the :ref:`doc.KDP`,
or alternatively will be installed if you follow the :ref:`doc.standalone_installation`, except for the ``kimvv`` package, which can be installed using
``pip install kimvv``.

If you have not requested your new properties to be added to OpenKIM, remember that ``kim-tools`` will automatically look for property definitions in the ``local-props`` and ``local_props`` subdirectories of the current working directory. If you wish to put them somewhere else,
you can point the environment variable ``KIM_PROPERTY_PATH`` to their location. ``kim-tools`` will expand any globs, including recursive ``**``. Conversely, if the properties
you are using are already in OpenKIM, do not use local property definition files.

.. _doc.example_script:

Example Script for Running a Crystal Genome Test Driver
=======================================================

:ref:`sphx_glr_auto_examples_CrystalGenomeASEExample__TD_000000654321_000_run.py`

.. _doc.curated_tests:

Curated Set of Test Cases
=========================

Because your Test Driver will run on a wide variety of crystals and interatomic potentials,
it is important to test it on a diverse selection of both. Below is a curated table
of test cases that you can query for. These are arranged in increasing symmetry order
(from triclinic to cubic). Here are some important considerations when choosing test cases for your
Driver:

* Make sure your driver runs with both Simulator Models and Portable Models. If your Test Driver uses LAMMPS,
  make sure it runs with models that use both ``atom_style atomic`` and ``atom_style charge``. A good way to
  do this is to test with any Portable Model, and a Buckingham or ReaxFF Simulator Model. Note that both of
  these are idiosyncratic -- ReaxFF are very slow and not great at reaching a low tolerance in static minimizaton,
  while Buckingham are much faster, but are even worse at static relaxation.
* Try to test on a variety of symmetries. If your Test Driver only computes scalar properties, it is less important
  to test on every crystal family, while for tensor properties, you should test robustness across symmetries
  more thoroughly.
* Note the unstable example using ``EAM_Dynamo_AcklandMendelevSrolovitz_2004_FeP__MO_884343146310_006`` in the table.
  If your computational protocol is such that a phase transition is possible, you should test with this example.
  Because ``kim-tools`` does not allow phase transitions at this time, the symmetry check should raise an error.


You can explore more prototypes at
http://aflow.org/prototype-encyclopedia/, but it is not guaranteed that OpenKIM
will have results or a compatible interatomic potential
(https://openkim.org/browse/models/by-species).

All models must be istalled before being used. Note that in a script using multiple models, you will need to re-instantiate
your ``TestDriver`` class with a new model each time.

.. todo::

   Find a good monoclinic example to add

.. csv-table::
   :header-rows: 1
   :file: structure_table.csv
   :widths: 10, 10, 40, 40
   :delim: tab
