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
or alternatively will be installed if you follow the :ref:`doc.standalone_installation`.

``kim-tools`` will automatically look for property definitions in the ``local-props`` and ``local_props`` subdirectories of the current working directory. If you wish to put them somewhere else,
you can point the environment variable ``KIM_PROPERTY_PATH`` to their location. ``kim-tools`` will expand any globs, including recursive ``**``.

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
(from triclinic to cubic).

You can explore more prototypes at
http://aflow.org/prototype-encyclopedia/, but it is not guaranteed that OpenKIM
will have results or a compatible interatomic potential
(https://openkim.org/browse/models/by-species).

Every time you use a new model, you will need to install the model and re-instantiate
your ``TestDriver`` class.

.. todo::

   Rhombohedral crystals are currently not supported. Skip the AB_hR26_148_a2f_b2f structure below.

.. csv-table::
   :header-rows: 1
   :file: structure_table.csv
   :widths: 10, 10, 40, 40
   :delim: tab

Commands to install all required models in the KDP using ``kimitems``, or outside the KDP using ``kim-api-collections-management``:

.. code-block:: console

   kimitems install -D Sim_LAMMPS_Buckingham_FreitasSantosColaco_2015_SiCaOAl__SM_154093256665_000
   kimitems install -D EDIP_LAMMPS_Marks_2000_C__MO_374144505645_000
   kimitems install -D MEAM_LAMMPS_FernandezPascuet_2014_U__MO_399431830125_002
   kimitems install -D SNAP_LiHuChen_2018_NiMo__MO_468686727341_000
   kimitems install -D Sim_LAMMPS_Buckingham_MatsuiAkaogi_1991_TiO__SM_690504433912_000
   kimitems install -D MEAM_LAMMPS_JeongParkDo_2018_PdAl__MO_616482358807_002
   kimitems install -D Sim_LAMMPS_Vashishta_BroughtonMeliVashishta_1997_SiO__SM_422553794879_000
   kimitems install -D MEAM_LAMMPS_KoJimLee_2012_FeP__MO_179420363944_002
   kimitems install -D Sim_LAMMPS_BOP_MurdickZhouWadley_2006_GaAs__SM_104202807866_001
   kimitems install -D Sim_LAMMPS_ReaxFF_BrugnoliMiyataniAkaji_SiCeNaClHO_2023__SM_282799919035_000
   kimitems install -D EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005

.. code-block:: console

   kim-api-collections-management install user Sim_LAMMPS_Buckingham_FreitasSantosColaco_2015_SiCaOAl__SM_154093256665_000
   kim-api-collections-management install user EDIP_LAMMPS_Marks_2000_C__MO_374144505645_000
   kim-api-collections-management install user MEAM_LAMMPS_FernandezPascuet_2014_U__MO_399431830125_002
   kim-api-collections-management install user SNAP_LiHuChen_2018_NiMo__MO_468686727341_000
   kim-api-collections-management install user Sim_LAMMPS_Buckingham_MatsuiAkaogi_1991_TiO__SM_690504433912_000
   kim-api-collections-management install user MEAM_LAMMPS_JeongParkDo_2018_PdAl__MO_616482358807_002
   kim-api-collections-management install user Sim_LAMMPS_Vashishta_BroughtonMeliVashishta_1997_SiO__SM_422553794879_000
   kim-api-collections-management install user MEAM_LAMMPS_KoJimLee_2012_FeP__MO_179420363944_002
   kim-api-collections-management install user Sim_LAMMPS_BOP_MurdickZhouWadley_2006_GaAs__SM_104202807866_001
   kim-api-collections-management install user Sim_LAMMPS_ReaxFF_BrugnoliMiyataniAkaji_SiCeNaClHO_2023__SM_282799919035_000
   kim-api-collections-management install user EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005
