================================
kimvv Python package integration
================================

.. contents:: Table of Contents

The final step in creating your Test Driver is to make sure it will be able to cleanly integrate into the
`kimvv <https://github.com/openkim/kimvv>`_ Python package for the materials science community to run on their
own outside of the OpenKIM pipeline. ``kimvv`` is built using `setuptools <https://setuptools.pypa.io>`_ and
published on PyPI at https://pypi.org/project/kimvv/. This page will walk you through the process to test
that your Test Driver will operate correctly independent of your local directory structure or environment
configuration.

The testing workflow will run your Test Driver on FCC gold (prototype label ``A_cF4_225_a``) using 3 models:

  * a KIM Portable Model: ``EAM_Dynamo_Ackland_1987_Au__MO_754413982908_001``
  * a KIM Simulator Model: ``Sim_LAMMPS_ADP_StarikovGordeevLysogorskiy_2020_SiAuAl__SM_113843830602_000``
  * a non-KIM ASE Calculator: ``ase.calculators.lj.LennardJones(sigma=2.42324, epsilon=2.30580)``

Make sure your Test Driver is able to run these simulations in the first place before beginning this process (except
for the non-KIM calculator if your Test Driver uses LAMMPS).


Setting up the packaging files
==============================

Next, you need to make sure your Test Driver is ready for packaging. Here are the files you may need to provide. They should be in the same directory
as your ``kimspec.edn`` file, considered the root directory of your Test Driver (one level above the ``test_driver`` directory).
Examples of all of them exist in the example Test Driver: |example_url|

    * ``README.rst``: Your Test Driver should include documentation. In the future, these will be automatically compiled into documentation for the
      `kimvv <https://github.com/openkim/kimvv>`_ package. We recommend using `.rst format <https://docutils.sourceforge.io/rst.html>`_ to make this task
      easier. Think about what you would like users to know about the usage of your Test Driver.
    * ``requirements.txt``: This is the requirements file for Python dependencies. It should include ``kim-tools`` at the minimum, as well as any other packages
      your Python code imports directly. Note that this does NOT include Python Standard Library imports (if you can find it on
      https://docs.python.org/3/library/index.html, it does not need to be in ``requirements.txt``). It will be included in the ``kimvv`` project as
      `dynamic metadata <https://setuptools.pypa.io/en/latest/userguide/pyproject_config.html#dynamic-metadata>`_, which supports a subset of the
      `PyPI Requirements File Format <https://pip.pypa.io/en/stable/reference/requirements-file-format/>`_ (``-c/-r/-e`` and other flags are not supported).
    * ``MANIFEST.in``: if your Test Driver requires non-Python files to operate (e.g. data files, pre-made images, etc.), you must declare them in this file
      for them to be included in the distribution. See the example Test Driver for a simple sample. Full documentation regarding ``MANIFEST.in`` is
      `here <https://setuptools.pypa.io/en/latest/userguide/miscellaneous.html>`_. The paths should be specified relative to the ``MANIFEST.in`` file,
      the ``kimvv`` script ``pre_setup.py`` will automatically edit the paths for incorporation into the Python package.


Running the tests
=================

A Github workflow for automatically testing that your Test Driver will work in the ``kimvv`` package is included in the example Test Driver:
https://github.com/openkim-hackathons/CrystalGenomeASEExample__TD_000000654321_000/blob/main/.github/workflows/test_kimvv.yml

If you based your Test Driver off this example, you will already have this file in your repo. You may need to click on "Actions" in your repo and
allow actions to run if the action has not been running already. If you created your repo from scratch, copy this workflow into your repo. It is set
up to run every time you push to the repo, or it can be triggered manually.

In the workflow file, you should edit the environment variable ``CRYSTAL_GENOME_TEST_ARGS``, which is a string literal containing a dictionary
of the additional arguments you would like to use for testing. The ASE example test driver has a simple example. A more elaborate example can be
found in the LAMMPS example Test Driver: https://github.com/openkim-hackathons/NPTCrystalStructure/blob/main/.github/workflows/test_kimvv.yml. Here
you can see various arguments set to make the molecular dynamics run in a reasonable time.

Your Test Driver must pass all the jobs in the workflow. Additionally, you should inspect the job output to see whether there is a
``WARNING Your Test Driver is unable to run with non-KIM calculators...`` during the "build and test kimvv" step. This warning is only acceptable
if you strictly need a KIM model name to run (e.g. if you are using LAMMPS). Otherwise, there should be no such warning.

.. todo::

  Consider adding more detail and pictures to this section if needed.
