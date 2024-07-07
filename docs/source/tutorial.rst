

=============================
Test Driver Creation Tutorial
=============================
  
This is a tutorial for creating OpenKIM Test Drivers in the Crystal Genome framework for arbitrary crystals. It follows the example Test Driver hosted at |example_url|. You are encouraged to use it as a template when writing your own Test Driver. 

.. note::
    If you are working in the KIM Developer Platform and intend to test the integration of your Test Driver into the OpenKIM pipeline, you should put your Test Driver in ``/home/openkim/test-drivers/<Extended KIM ID Prefix>__TD_DDDDDDDDDDDDD_VVV``, where ``<Extended KIM ID Prefix>`` is a short alphanumeric description of your test driver, ``DDDDDDDDDDDD`` is an arbitrary unique 12-digit integer code (it is conventional for in-development codes to start with `000000`), and ``VVV`` is the version starting with ``000``. For more info see https://openkim.org/doc/schema/kim-ids/.

.. toctree::
   :maxdepth: 2
   :caption: Steps for Creating a Crystal Genome Test Driver

   tutorial_property
   test_driver
   auto_examples/CrystalGenomeASEExample__TD_000000654321_000/debug
   tutorial_pipeline
