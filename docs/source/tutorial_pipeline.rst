================================================
Metadata and KIM Processing Pipeline Integration
================================================

.. contents:: Table of Contents

The final step in creating your Test Driver is to set up the auxiliary files needed for the KIM Processing Pipeline to use it.
As with the rest of this documentation, this page follows the example test driver from |example_url|. Most of the files are
pre-configured to work as-is for most Test Drivers, and only a few changes are needed.

To test that everything is correctly configured, you will need to use the :ref:`doc.KDP` and place your Driver in the 
:ref:`correct path <note.td_path>`.

.. todo::
    
    In the future, we will automate the manual processes described here.

First, let's go over the files you will need to change.

``test_generator.json``
=======================

This is the file that specifies the inputs to your Test Driver. It contains all structures currently in Crystal Genome.
Here are the first few lines of it. Each line is a JSON dictionary. Each line corresponds to a separate `KIM Test <https://openkim.org/doc/evaluation/kim-tests/>`_
that will be created using your Test Driver.

.. literalinclude:: ../../examples/CrystalGenomeASEExample__TD_000000654321_000/test_generator.json
    :language: json
    :lines: 1-4

As you can see, it specifies the crystal structure and allows you to pass any additional arguments you defined in your ``_calculate()`` method as a sub-dictionary
under ``crystal_genome_test_args``. At the very least, you will need to edit this sub-dictionary to reflect the inputs your Test Driver uses. You can omit any
arguments you wish to set to their default value, and if your ``_calculate()`` takes no additional arguments, you may make it an empty dictionary: ``{}``. If your
Test Driver does not use temperature and/or stress, you may omit these keys as well (the example Test Driver does not use them, but they are left for demonstration).

If you wish to create Tests for different inputs (e.g. multiple temperatures), you must create a separate line for each one. For example, if we wanted to create Tests
using the example Test Driver that set ``max_volume_scale=0.1`` and different Tests that set ``max_volume_scale=0.2``, both for all structures, we would double the number
of lines in ``test_generator.json``.

.. todo::
    
    Currently the temperature and stress inputs only set the internal values of those variables in the :class:`~kim_test_utils.test_driver.CrystalGenomeTestDriver` class, and do not affect querying. The structure queried for is always the zero temperature, zero pressure structure

Generating and Running Tests
============================

There is one more file that you will need to change from the example Test Driver, but it will be easier to contextualize if we first demonstrate how Tests are generated and run in the KIM Developer Platform using the example Test Driver.

To generate tests from the example Test Driver you've placed in ``~/test-drivers/CrystalGenomeASEExample__TD_000000654321_000`` of your KDP, run the following command. Note that it takes several minutes to generate Tests for the large number of lines in ``test_generator.json``! You can either stop the process (Ctrl-C or Command-C) after a few tests have been generated, or edit ``test_generator.json`` to only leave a few lines remaining.

.. code-block:: bash

    kimgenie tests --add-random-kimnums --test-driver CrystalGenomeASEExample__TD_000000654321_000

.. note:: 
    The ``--add-random-kimnums`` option has added a ``kimnum`` key to each line in ``test_generator.json``. Do not use this option after these keys have been added.

Assuming you have left the first few lines of ``test_generator.json`` intact, you will have created a test for FCC silver. Let's install a model for silver and run it with the created test (adding the custom property as well if we have not done so already). The random ``kimnum`` generated for your test will vary, but the ``pipeline-run-pair`` command accepts wildcards:

.. code-block:: bash

    kimitems install -D  EAM_Dynamo_AcklandTichyVitek_1987_Ag__MO_212700056563_005
    add_or_update_property ~/test-drivers/CrystalGenomeASEExample__TD_000000654321_000/local-props/energy-vs-volume-isotropic-crystal.edn
    pipeline-run-pair pipeline-run-pair CrystalGenomeASEExample_A_cF4_225_a_Ag_0_1__TE_* EAM_Dynamo_AcklandTichyVitek_1987_Ag__MO_212700056563_005 -v

You should see the output of the test, and there should be a new directory in ``~/test-results/`` with the ``results.edn`` file containing the resulting KIM Property Instance, 
as well as several other files documenting the run.

For more info regarding the utilities available in the KIM Developer Platform, see https://openkim.org/doc/evaluation/kim-developer-platform/, consult the ``README.txt`` in the home directory of the KDP, or run any of them with the ``-h`` option.

``test_template/kimspec.edn.genie``
===================================

If you look in your ``~/tests/`` directory, you will see that a subdirectory has been created for each Test. Within it is a file named ``kimspec.edn``, which is a metadata file required to be provided with every piece of KIM content. 
The ``kimgenie`` utility generates these from the ``test_template/kimspec.edn.genie`` file in the Test Driver directory using Jinja2 templating. 
A Jinja2 tutorial is available here: https://ttl255.com/jinja2-tutorial-part-1-introduction-and-variable-substitution/.
Below is the file in from the example Test Driver. Compare with one of the rendered ``kimspec.edn`` files in ``~/tests/`` to see how the templating works.

.. literalinclude:: ../../examples/CrystalGenomeASEExample__TD_000000654321_000/test_template/kimspec.edn.genie
    :language: jinja

Here are the keys you will need to edit. You can leave the templating involving ``TEST_DRIVER_NAME``, ``prototype_label``, ``stoichiometric_species`` and ``kimnum`` as-is, 
but you will need to change descriptive prose, as well as render the arguments specific to your Test Driver (if any). Note, for example, how the ``max_volume_scale`` argument is templated
in the fields below. See https://openkim.org/doc/schema/kimspec/ for a detailed definition of each field.

    * ``title`` and ``description``: change these to be relevant to your Test Driver and descriptive of the individual Tests.
    * ``extended-id``: See https://openkim.org/doc/schema/kim-ids/ for recommendations. Note how a ``replace`` operation is used on ``crystal_genome_test_args.max_volume_scale``, as ``.`` is not an allowed character.
    * ``developer``: List of all developers' OpenKIM User IDs. To get a User ID, sign up for an account at https://openkim.org/new-account. You may also add an additional optional key ``implementer`` for anyone who contributed only programming, not core intellectial content. The two lists should not share any entries.

In most cases, this should be all you need to change. Use ``kimgenie`` to generate some Tests and check that ``kimspec.edn`` looks correct.

Other files
===========

These are the other required files in a Test Driver. In most cases, other than the README, you can copy them from the example Test Driver as-is.

    * ``README.txt``: Your Test Driver should include some kind of documentation, but the format is up to you.
    * ``kimspec.edn``: This is the metadata for the Test Driver itself. You may leave this in the minimal form provided in the example Test Driver, or you can fill in some of the fields following https://openkim.org/doc/schema/kimspec/. Any missing fields will be populated by the Web form when you :ref:`submit your Driver <doc.submit>`.
    * ``Makefile``: Unless your Test Driver requires compilation, leave this dummy Makefile untouched.
    * ``runner``: This is the executable that invokes the ``TestDriver`` class when running in the KIM Pipeline.
    * ``test_template/runner``: This is a wrapper executable in each Test that invokes the Test Driver's ``runner``. You should never have to change this file.
    * ``test_template/pipeline.stdin.tpl.genie``: This is the Jinja2 template file for passing inputs to the ``runner``.
    * ``test_template/dependencies.edn.genie``: This Jinja2 template file specifies the Tests that the generated Tests are dependent on.
    * ``Makefile``, ``LICENCE``: Always leave these as-is.
