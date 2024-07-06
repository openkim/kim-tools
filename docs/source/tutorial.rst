=============================
Test Driver Creation Tutorial
=============================
  
This is a tutorial for creating OpenKIM Test Drivers in the Crystal Genome framework for arbitrary crystals. It follows the example Test Driver hosted at https://github.com/openkim-hackathons/CrystalGenomeASEExample__TD_000000654321_000. You are encouraged to use it as a template when writing your own Test Driver. 

.. note::
    If you are working in the KIM Developer Platform and intend to test the integration of your Test Driver into the OpenKIM pipeline, you should put your Test Driver in ``/home/openkim/test-drivers/<Extended KIM ID Prefix>__TD_DDDDDDDDDDDDD_VVV``, where ``<Extended KIM ID Prefix>`` is a short alphanumeric description of your test driver, ``DDDDDDDDDDDD`` is an arbitrary unique 12-digit integer code (it is conventional for in-development codes to start with `000000`), and ``VVV`` is the version starting with ``000``. For more info see https://openkim.org/doc/schema/kim-ids/.

.. _doc.tutorial.property:

Property Definition
===================

.. toctree::
   :maxdepth: 2
     
   tutorial_property

The first step is to create an OpenKIM *Property Definition* to describe the material property your Test Driver will calculate. 

Learn how to do so here: :doc:`tutorial_property`. Once you have created your property definition file, you need to run the ``add_or_update_this_property`` command-line utility included with this package. Pass the path to the property definition file as a command line argument. For example, to add the property provided with the ``CrystalGenomeASEExample__TD_000000654321_000``, assuming you have placed the driver into ``/home/openkim/test-drivers/``, run

.. code-block:: bash

    add_or_update_this_property /home/openkim/test-drivers/CrystalGenomeASEExample__TD_000000654321_000/local-props/energy-vs-volume-isotropic-crystal.edn

Main Computational Code
=======================

.. toctree::
   :maxdepth: 2
     
   auto_examples/CrystalGenomeASEExample__TD_000000654321_000/test_driver

The next step is to write the main computational code. It must be contained in a Python file named ``test_driver.py``, although you may include as many other Python files as you wish for utility functions. The file included in the example Test Driver doubles as a tutorial and can be found here: :doc:`auto_examples/CrystalGenomeASEExample__TD_000000654321_000/test_driver`

Metadata and KIM Processing Pipeline Integration
================================================

.. toctree::
    :maxdepth: 2

    tutorial_pipeline


