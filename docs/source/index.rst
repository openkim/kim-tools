.. kim-tools documentation master file, created by
   sphinx-quickstart on Mon Jun 17 01:51:39 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==========================================
Welcome to kim-tools's documentation!
==========================================

.. figure:: underconstruction.gif

   :sub:`image source: http://textfiles.com/underconstruction/`

.. toctree::
   :maxdepth: 3
   :caption: Contents:
   
   tutorial
   modules

Introduction
============

This package is designed to facilitate writing `OpenKIM Test Drivers <https://openkim.org/doc/evaluation/kim-tests/>`_ in the Crystal Genome framework for arbitrary crystals. Fundamentally, a KIM Test Driver is an arbitrary executable named ``runner`` that takes custom inputs from ``stdin`` using custom and standard Jinja2 templating, and writes a file named ``output/results.edn`` containing one or more `KIM Property Instances <https://openkim.org/doc/schema/properties-framework/>`_. This minimal specification is all the *KIM Processing Pipeline* needs to run KIM Tests and allows a great amount of flexibility for developers. However, the open-endedness makes creation of new Test Drivers difficult to tutorialize, and the resulting Test Drivers have a non-uniform input syntax, making them cumbersome to use outside of the KIM Pipeline.

*Crystal Genome* is the OpenKIM thrust for generalizing our testing framework to apply to arbitrary crystals. It uses AFLOW Prototype Labels (see figure below) to describe arbitrary crystal structures in a complete, concise, and human-readable manner. A thorough description of the AFLOW prototype label can be found in Part IIB here: https://arxiv.org/abs/2401.06875. Because it is unreasonable for every OpenKIM test developer to learn the crystallography required to write Crystal Genome tests from scratch, common input, output, and analysis tasks must be easy to use automatically.

.. figure:: aflow.png
   :align: center
   :width: 500px

   The AFLOW prototype label. In short, it is a string providing the stoichiometry, Pearson symbol, space group, and Wyckoff positions of an arbitrary crystal. Image source and more info: https://arxiv.org/pdf/2401.06875. 


This package addresses these issues by providing a base class :class:`~kim_tools.CrystalGenomeTestDriver` (itself inheriting from a yet more general class :class:`~kim_tools.KIMTestDriver`) that automates many of the common programming tasks needed to write an OpenKIM Test Driver for Crystal Genome. Additionally, due to the common interface, Test Drivers written using this base class are easier to invoke outside of the OpenKIM infrastructure, and can work with arbitrary ASE :class:`~ase.Atoms` objects. If a Test Driver uses only ASE for computations, it can even work with arbitrary ASE :class:`~ase.calculators.calculator.Calculator` objects. In the future, all ASE-only Test Drivers written using :mod:`kim_tools` will be incorporated into the `kimvv <https://github.com/openkim/kimvv>`_ package to be distributed for users to test their own models.

Contact us at https://openkim.org/contact/ with any comments or questions.

Test Driver Creation Tutorial
=============================

If you wish to use this package to develop an OpenKIM Test Driver, after setting up a suitable development environment (detailed on the rest of this page), follow the :doc:`tutorial`. 

.. _doc.KDP:

KIM Developer Platform
======================

The `KIM Developer Platform (KDP) <https://openkim.org/doc/evaluation/kim-developer-platform/>`_ is a Docker image providing an emulation of the KIM Processing Pipeline. It is the recommended environment for developing OpenKIM content, and is required to fully test integration of Test Drivers into the pipeline before they are submitted to `openkim.org <https://openkim.org/>`_. ``kim-tools`` and all requirements are included in KDP version 1.4.0 and higher.

.. _doc.standalone_installation:

Standalone Installation
=======================

Standalone usage of ``kim-tools`` is possible and makes sense if you wish to use the package on an HPC resource, as most do not support running Docker images. 

-------------
Prerequisites
-------------
Before installing the package using ``pip``, you will need to install the KIM API and the AFLOW software. It is likely that you will have to do some configuration for your specific machine for things like Conda instalaltion directories, ``pip`` installation directories, etc.

KIM API
-------
Installation instructions for the KIM API can be found here: https://openkim.org/doc/usage/obtaining-models/. Recommended installation methods are building from source or installing from conda-forge. If you need `LAMMPS <https://www.lammps.org>`_ as well (because you are running or developing Test Drivers that use it, or you wish to use LAMMPS `Simulator Models <https://openkim.org/doc/repository/kim-content/>`_), a good option is to install the ``lammps`` package from conda-forge which includes ``kim-api`` as a dependency.

There may be issues mixing ``pip`` and Conda. I did not encounter any arising directly from ``kim-tools``, but when developing a test that uses ``numdifftools``, I found that it was only possible to install through Conda, not ``pip``, for example.

AFLOW
-----
The installation script for ``kim-tools`` will check that the ``aflow`` executable is in your ``PATH`` and working correctly. Installation instructions are here: http://aflow.org/install-aflow/. If the automatic installers/downloaders do not work for you, I have found building from source manually to be straightforward and well-behaved. Simply download the source code, issue ``make``, and copy the resulting executable into your ``PATH`` manually or using ``make install``. An executable file named ``aflow_data`` must also be present in your ``PATH``. For the purposes of ``kim-tools``, this may be an empty file with execute permissions, it does not need to be the actual ``aflow_data`` file built from the AFLOW source.

Installation
------------
Once the above prerequisites are satisfiied, you can install using ``pip``:

.. code-block:: bash

   pip install --user git+https://github.com/openkim/kim-tools.git

If you are installing on an HPC cluster, you may need to use different options to achieve the desired result.

API Documentation, Indices and tables
=====================================

* :ref:`modindex`
* :ref:`search`
* :ref:`genindex`
