=====================================================
Testing and Debugging Your Crystal Genome Test Driver
=====================================================

When your Test Driver is integrated into the OpenKIM Processing Pipeline, the auxiliary files supplied with your
Test Driver will be used to pass inputs and parse outputs from the calculation you have implemented in ``test_driver.py``.
This functionality is also available for testing in the KIM Developer Platform. More information about this is avalable
here: :doc:`../../tutorial_pipeline`, however it is likely that you will first wish to test your Driver by invoking it directly.

Examples for how to do that are found in the script ``run.py`` from |example_url|, documented below in the :ref:`example script <sphx_glr_auto_examples_CrystalGenomeASEExample__TD_000000654321_000_run.py>` below.
A practical guide for testing your Driver on a variety of crystal structures is found at the bottom of the page in the section ':ref:`doc.curated_tests`'.
Everything you need to run the example script is containerized in the :ref:`doc.KDP`,
or alternatively will be installed if you follow the :ref:`doc.standalone_installation`. 

``kim-tools`` will automatically look for property definitions in the ``local-props`` and ``local_props`` subdirectories of the current working directory. If you wish to put them somewhere else,
you can point the environment variable ``KIM_PROPERTY_PATH`` to their location. ``kim-tools`` will expand any globs, including recursive ``**``.

Example Script for Running a Crystal Genome Test Driver
=======================================================

:ref:`sphx_glr_auto_examples_CrystalGenomeASEExample__TD_000000654321_000_run.py`

.. _doc.curated_tests:

Curated Set of Test Cases
=========================

Because your Test Driver will run on a wide variety of crystals and interatomic potentials, 
it is important to test it on a diverse selection of both. Below is a curated table
of test cases that you can query for. These are arranged in increasing symmetry order
(from triclinic to cubic).

You can explore more prototypes at 
http://aflow.org/prototype-encyclopedia/, but it is not guaranteed that OpenKIM
will have results or a compatible interatomic potential 
(https://openkim.org/browse/models/by-species).

Every time you use a new model, you will need to install the model and re-instantiate
your ``TestDriver`` class.

.. csv-table:: 
   :header-rows: 1
   :file: structure_table.csv
   :widths: 10, 10, 40, 40
   :delim: tab