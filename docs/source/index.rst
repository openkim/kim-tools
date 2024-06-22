.. kim-test-utils documentation master file, created by
   sphinx-quickstart on Mon Jun 17 01:51:39 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to kim-test-utils's documentation!
==========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Introduction
------------

This package is designed to facilitate writing `OpenKIM Test Drivers <https://openkim.org/doc/evaluation/kim-tests/>`_, especially Test Drivers in the Crystal Genome framework for arbitrary crystals. Fundamentally, a KIM Test Driver is an arbitrary executable named ``runner`` that takes custom inputs from ``stdin`` using custom and standard Jinja2 templating, and writes a file named ``output/results.edn`` containing one or more `KIM Property Instances <https://openkim.org/doc/schema/properties-framework/>`_. This minimal specification is all the *KIM Processing Pipeline* needs to run KIM Tests and allows a great amount of flexibility for developers. However, the open-endedness makes creation of new Test Drivers difficult to tutorialize, and the resulting Test Drivers have a non-uniform input syntax, making them cumbersome to use outside of the KIM Pipeline.

This package addresses these issues by providing a base class :class:`CrystalGenomeTestDriver`

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
