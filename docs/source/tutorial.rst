

=============================
Test Driver Creation Tutorial
=============================
  
This is a tutorial for creating OpenKIM Test Drivers in the Crystal Genome framework for arbitrary crystals. It follows the example Test Driver hosted at |example_url|. You are encouraged to use it as a template when writing your own Test Driver. For example, if you are working in the :ref:`doc.KDP`, you can clone the example Test Driver with the following shell commands:

.. parsed-literal::

    cd ~/test-drivers
    git clone \ |example_url|\ 

.. _note.td_path:
.. note::
    If you are working in the KIM Developer Platform and intend to test the integration of your Test Driver into the OpenKIM pipeline, the location and directory name of the Test Driver is important.
    This means that when you should follow the above instructions exactly to clone into the correct location, and when you create your own Test Driver, you should place it in ``/home/openkim/test-drivers/<Extended KIM ID Prefix>__TD_DDDDDDDDDDDDD_VVV``. ``<Extended KIM ID Prefix>`` is a short alphanumeric description of your test driver, ``DDDDDDDDDDDD`` is an arbitrary unique 12-digit integer code (it is conventional for in-development codes to start with `000000`), and ``VVV`` is the version starting with ``000``. For more info see https://openkim.org/doc/schema/kim-ids/.

.. toctree::
   :maxdepth: 2
   :caption: Steps for Creating a Crystal Genome Test Driver

   tutorial_property
   test_driver
   auto_examples/CrystalGenomeASEExample__TD_000000654321_000/debug
   tutorial_pipeline

.. _doc.submit:

Submitting Your Test Driver
===========================

When your Test Driver is complete and fully tested, take the following final steps to get your Test Driver submitted to OpenKIM.org and running in the Pipeline:

    * Request a member of the OpenKIM team to permanently add your Property Definition to the collection of `KIM Property Definitions <https://openkim.org/properties>`_.
    * Create an OpenKIM account at https://openkim.org/new-account.
    * Submit a tarball of your Test Driver at https://openkim.org/contribute/test-driver/. We recommend issuing a release of your repository using Git or Github, the ``.gitattributes`` file provided with the example Test Driver already has the correct settings to only export the needed files. Follow the instructions in the form to fill out the information about your Test Driver.
