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

  * a KIM Portable Model: ``LennardJones612_UniversalShifted__MO_959249795837_003``
  * a KIM Simulator Model: ``Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu``
  * a non-KIM ASE Calculator: ``ase.calculators.lj.LennardJones(sigma=2.42324, epsilon=2.30580, rc=9.69298)``

Make sure your Test Driver is able to run these simulations in the first place before beginning this process (except
for the non-KIM calculator if your Test Driver uses LAMMPS). The two KIM models are example models included in all
builds of the KIM API by default, except the conda-forge distribution for Apple Silicon Mac.

Note that because these models are example models that are not on OpenKIM.org, querying for their equilibrium structures
will not work. Instead you should start with the relaxed structure using ``kimvv.EquilibriumCrystalStructure``. This is what
the testing workflow will do:

.. code-block:: python

  from ase.build import bulk
  from kimvv import EquilibriumCrystalStructure

  atoms_init = bulk('Au')

  # Instantiate the Equilibrium Driver with your model
  kim_model_name = "Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu"
  ecs = EquilibriumCrystalStructure(kim_model_name)

  # Relax the structure. ECS will return multiple properties, any of them will do as they all contain the
  # crystal description
  relaxed_structure = ecs(atoms_init)[0]

  # Run your TD with `relaxed_structure` as the input


Forking the ``kimvv`` repo
==========================

.. list-table::
   :class: borderless

   * - First, make a fork of the ``kimvv`` repo under your own Github account by navigating to https://github.com/openkim/kimvv and clicking the "fork" button:
     - .. image:: fork_kimvv.png
        :width: 500px
   * - Once the fork is in your account, before making any changes, allow Github actions to run by clicking on the "actions" button and agreeing to run actions
       on the resulting page:
     - .. image:: actions.png
        :width: 500px

Setting up the packaging files
==============================

Next, you need to make sure your Test Driver is ready for packaging. Here are the files you may need to provide. They should be in the same directory
as your ``kimspec.edn`` file, considered the root directory of your Test Driver (one level above the ``test_driver`` directory).
Examples of all of them exist in the example Test Driver: |example_url|

    * ``README.rst``: Your Test Driver should include documentation. In the future, these will be automatically compiled into documentation for the
      `kimvv <https://github.com/openkim/kimvv>`_ package. We recommend using `.rst format <https://docutils.sourceforge.io/rst.html>`_ to make this task
      easier. Think about what you would like users to know about the usage of your Test Driver.
    * ``requirements.txt``: This is the requirements file for Python dependencies. It should include ``kim-tools`` at the minimum, as well as any other packages
      your Python code imports directly. It will be included in the ``kimvv`` project as
      `dynamic metadata <https://setuptools.pypa.io/en/latest/userguide/pyproject_config.html#dynamic-metadata>`_, which supports a subset of the
      `PyPI Requirements File Format <https://pip.pypa.io/en/stable/reference/requirements-file-format/>`_ (``-c/-r/-e`` and other flags are not supported).
    * ``MANIFEST.in``: if your Test Driver requires non-Python files to operate (e.g. data files, pre-made images, etc.), you must declare them in this file
      for them to be included in the distribution. See the example Test Driver for a simple sample. Full documentation regarding ``MANIFEST.in`` is
      `here <https://setuptools.pypa.io/en/latest/userguide/miscellaneous.html>`_. The paths should be specified relative to the ``MANIFEST.in`` file,
      the ``kimvv`` script ``pre_setup.py`` will automatically edit the paths for incorporation into the Python package.


Issuing a release and adding it to your fork of ``kimvv``
=========================================================

Next, you need to issue a release of your Test Driver. Because the release that is fully tested and ready for OpenKIM submission
will be named and tagged "v000", it is recommended that you tag your release something like "v000b0". Make sure your ``.gitattributes``
file is correctly configured to ignore any files and directories that should be excluded from releases
(e.g. ``run.py``, ``output/``, ``local-props/``). See the
`example .gitattributes <https://github.com/openkim-hackathons/CrystalGenomeASEExample__TD_000000654321_000/blob/794664404260f9a6fc556e9401dba4851cdeb9c5/.gitattributes>`_.

.. list-table::
   :class: borderless

   * - Issue a release by clicking the "Create a new release" link in your Test Driver repository:
     - .. image:: release_start.png
        :width: 500px
   * - Once the release is issued, you will need the URL of the ``.tar.gz`` file created by GitHub as part of it:
     - .. image:: release_download.png
        :width: 500px

Next, you need to add the URL of the ``.tar.gz`` to the ``pre_setup.py`` script in your fork of ``kimvv``. An example fork for testing the example Test Driver is available,
and here is where you add the URL: `pre_setup.py <https://github.com/openkim-hackathons/kimvv-example-driver-testing-fork/blob/main/pre_setup.py>`_.
Finally, you need to add your Property Definitions to the ``test/local-props`` directory of your ``kimvv`` fork. Alternatively, you can skip this by requesting an OpenKIM team member
to already publish your Property and add it to the ``kim-property`` package.

Checking that your fork passes the tests
========================================

After making the above changes, commit and push the changes to your fork of ``kimvv``. The Github Actions workflow will automatically check that the ``kimvv`` package
containing your driver correctly installs and runs the previously described simulations.

.. list-table::
   :class: borderless

   * - When you click on "Actions" in your ``kimvv`` fork on Github after pushing your changes, you should see the test workflow in progress, failed, or succeeded.
       the example on the right shows a failed and successful run:
     - .. image:: actions_results.png
        :width: 500px
   * - If your Test Driver uses LAMMPS, a successful run is all that is required. If it does not, you should take one additional step to see that your Test Driver
       runs with non-KIM ASE calculators. Click on the successful run, and click on one of the jobs (e.g. ``test (ubuntu-22.04")``). Expand the "Run tests" section
       and check that there is no ``WARNING Your Test Driver is unable to run with non-KIM calculators...``.
     - .. image:: actions_warnings.png
        :width: 500px
