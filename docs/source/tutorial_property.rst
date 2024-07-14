==================================
Creating a Crystal Genome Property
==================================

The first step is to create an OpenKIM *Property Definition* to describe the material property your Test Driver will calculate. 

Under the `KIM Properties Framework <https://openkim.org/doc/schema/properties-framework/>`_, a *Property Definition* is a file in `EDN format <https://openkim.org/doc/schema/edn-format/>`_ that defines the fields required to fully characterize a material property.

An example Property Definition is provided with the example Test Driver |example_url| at ``local-props/energy-vs-volume-isotropic-crystal.edn`` and is shown at the bottom of this page. All properties in the Crystal Genome framework share a set of keys for describing the crystal using the AFLOW Prototype Label, as well as other metadata. See the comments in the example Property Definition to see which sections you should change and which you should retain.

We strive to make all Crystal Genome properties as universally applicable to arbitrary crystals as possible. You should think carefully about what is needed for a minimal, yet complete, description of your material property when applied to a generic crystal with arbitrary symmetry. To demonstrate this point, consider that the example property here is *NOT* suitable for general crystals. This is because isotropic expansion and contraction is an arbitrary displacement boundary condition for any non-cubic crystal not corresponding to any specific loading. For cubic crystals it is a meaningful property, as it always corresponds to hydrostatic stress -- see the non-Crystal Genome version of the property here: https://openkim.org/properties/show/2014-04-15/staff@noreply.openkim.org/cohesive-energy-relation-cubic-crystal. You are encouraged to work with the OpenKIM team to develop your property definition (contact us individually or at https://openkim.org/contact/).

Here are some examples of existing Crystal Genome properties:
    * https://openkim.org/properties/show/2024-07-10/staff@noreply.openkim.org/bulk-modulus-isothermal-npt
    * https://openkim.org/properties/show/2024-07-10/staff@noreply.openkim.org/elastic-constants-isothermal-npt

It is possible that your Test Driver will write multiple material properties. In this case, create a separate property definition file for each. For example, the bulk modulus and elastic constants shown above are separate properties.

Once you have created your property definition file, you need to run the ``add_or_update_this_property`` command-line utility included with this package. Pass the path to the property definition file as a command line argument. For example, to add the property provided with ``CrystalGenomeASEExample__TD_000000654321_000``, assuming you have placed the driver into ``/home/openkim/test-drivers/``, run 

.. code-block:: bash

    add_or_update_this_property /home/openkim/test-drivers/CrystalGenomeASEExample__TD_000000654321_000/local-props/energy-vs-volume-isotropic-crystal.edn

Once your property is finalized, you will need to request a member of the OpenKIM team to permanently add it to the collection of `KIM Property Definitions <https://openkim.org/properties>`_. Your property should only be added after you have gone through the rest of the steps in this tutorial to make sure you do not encounter any unanticipated problems while developing your Test Driver. However, you are encouraged to contact us early in the development to begin iterating on your property definiton.

Tips
====

* The `KIM Properties Framework <https://openkim.org/doc/schema/properties-framework/>`_ allows for flexible (optional) reporting of uncertainty information for every key, so you do not need to add uncertainties as separate keys. See how to do so at the end of the :ref:`doc.example_test_driver`.
* Do not add pressure as a key for constant-pressure simulations, report it as a hydrostatic stress state in the ``cell-cauchy-stress`` key.

Example Property
================

.. literalinclude:: ../../examples/CrystalGenomeASEExample__TD_000000654321_000/local-props/energy-vs-volume-isotropic-crystal.edn
    :language: clojure