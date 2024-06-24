.. kim-test-utils documentation master file, created by
   sphinx-quickstart on Mon Jun 17 01:51:39 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==========================================
Welcome to kim-test-utils's documentation!
==========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. figure:: underconstruction.gif

   :sub:`image source: http://textfiles.com/underconstruction/`

Introduction
============

This package is designed to facilitate writing `OpenKIM Test Drivers <https://openkim.org/doc/evaluation/kim-tests/>`_ in the Crystal Genome framework for arbitrary crystals. Fundamentally, a KIM Test Driver is an arbitrary executable named ``runner`` that takes custom inputs from ``stdin`` using custom and standard Jinja2 templating, and writes a file named ``output/results.edn`` containing one or more `KIM Property Instances <https://openkim.org/doc/schema/properties-framework/>`_. This minimal specification is all the *KIM Processing Pipeline* needs to run KIM Tests and allows a great amount of flexibility for developers. However, the open-endedness makes creation of new Test Drivers difficult to tutorialize, and the resulting Test Drivers have a non-uniform input syntax, making them cumbersome to use outside of the KIM Pipeline.

This package addresses these issues by providing a base class :class:`~kim_test_utils.CrystalGenomeTestDriver` that automates many of the common programming tasks needed to write an OpenKIM Test Driver. Additionally, due to the common interface, tests written using this base class are easier to invoke outside of the OpenKIM infrastructure, and in fact should work with arbitrary ASE :class:`~ase.Atoms` and :class:`~ase.calculators.calculator.Calculator` objects. In the future, all tests written using this package will be incorporated into the ``kimvv`` package to be distributed for users to test their own models.

KIM Developer Platform
======================

The `KIM Developer Platform (KDP) <https://openkim.org/doc/evaluation/kim-developer-platform/>`_ is a Docker image providing an emulation of the KIM Processing Pipeline. It is the recommended environment for developing OpenKIM content, and is required to fully test integration of Test Drivers into the pipeline before they are submitted to `openkim.org <https://openkim.org/>`_. ``kim-test-utils`` and all requirements are included in KDP version 1.4.0 and higher.

Standalone Installation
=======================

Standalone usage of ``kim-test-utils`` is possible and makes sense if you wish to use the package on an HPC resource, as most do not support running Docker images. You will likely need to troubleshoot the 

-------------
Prerequisites
-------------
Before installing the package using ``pip``, you will need to install the KIM API and the AFLOW software. It is likely that you will have to do some configuration for your specific machine for things like Conda instalaltion directories, ``pip`` installation directories, etc.

KIM API
-------
Installation instructions for the KIM API can be found here: https://openkim.org/doc/usage/obtaining-models/. Recommended installation methods are building from source or installing from conda-forge. If you need `LAMMPS <https://www.lammps.org>`_ as well (because you are running or developing Test Drivers that use it, or you wish to use LAMMPS `Simulator Models <https://openkim.org/doc/repository/kim-content/>`_), a good option is to install the ``lammps`` package from conda-forge which includes ``kim-api`` as a dependency.

There may be issues mixing ``pip`` and Conda. I did not encounter any arising directly from ``kim-test-utils``, but when developing a test that uses ``numdifftools``, I found that it was only possible to install through Conda, not ``pip``, for example.

AFLOW
-----
The installation script for ``kim-test-utils`` will check that the ``aflow`` executable is in your ``PATH`` and working correctly. Installation instructions are here: http://aflow.org/install-aflow/. If the automatic installers/downloaders do not work for you, I have found building from source manually to be straightforward and well-behaved. Simply download the source code, issue ``make``, and copy the resulting executable into your ``PATH`` manually or using ``make install``. An executable file named ``aflow_data`` must also be present in your ``PATH``. For the purposes of ``kim-test-utils``, this may be an empty file with execute permissions, it does not need to be the actual ``aflow_data`` file built from the AFLOW source.

Installation
------------
Once the above prerequisites are satisfiied, you can install using ``pip``:

.. code-block:: bash

   pip install --user git+https://github.com/openkim/kim-test-utils.git

If you are installing on an HPC cluster, you may need to use different options to achieve the desired result. Note that if you wish to be able to use the included command-line utility ``add_or_update_this_property`` (you will if you are developing a new Crystal Genome Test Driver), the prerequisite package ``kim-property`` must be installed in an editable location. The ``--user`` option included in the example command does so.

Crystal Genome Test Driver Tutorial
===================================

:doc:`auto_examples/index`

API Documentation, Indices and tables
=====================================

* Autogenerated API documentation is here: :ref:`modindex`
* :ref:`search`
* :ref:`genindex`
