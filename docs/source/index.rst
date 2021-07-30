.. Genonets documentation master file, created by
   sphinx-quickstart on Tue Jul 27 12:43:12 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Genonets' documentation!
===================================

This package provides a high level interface for construction and analysis of
genotype networks from data.

Following features are available:

* Parsing of genotype-phenotype maps from an input file, provided in the
  :doc:`genonets input-file format <input_format>`
* Creation of genotype networks from input data
* Various analyses on the constructed genotype networks
* Generation of result files with attributes from genotype-network-level
  analyses, as well as genotype-level analyses
* Creation of a phenotype network that shows evolvability and accessibility
  relationships between the genotype sets
* Generation of GML files corresponding to the created genotype networks and
  the phenotype network

.. _installation:

Installation
------------

Before installing Genonets, please use the following command to upgrade
:code:`pip` and :code:`setuptools`:

.. code-block:: console

   $ pip install -U pip setuptools

Then install Genonets using the following command:

.. code-block:: console

   $ pip install git+https://github.com/fkhalid/genonets.git

.. note::

   Only :code:`Python 3.8` and above are supported. Also, :code:`Linux` and
   :code:`MacOS` are supported; we plan to add support for :code:`Windows` at a
   later time.

Getting started
---------------

Tutorial I: The one-liner
^^^^^^^^^^^^^^^^^^^^^^^^^

The following code snippet shows the simplest possible way of using Genonets:

.. code-block:: python
   :linenos:

   from genonets.cmdl_handler import CmdParser
   from genonets.interface import Genonets

   Genonets(CmdParser().get_args(), process=True)

Yes, that's it. Two import statements, and just one line of code to create,
analyze, and save, all the genotype networks in the input file, as well the
phenotype network. In fact, these are the contents of the
:code:`minimal.py` sample file included in the package.

Assuming you have already installed Genonets and downloaded the sample input
file available :download:`here <sample_data/input.csv>`, the :code:`minimal.py`
program can be used to process the sample input file as follows.

.. code-block:: console

   $ python -u -m genonets.sample.minimal --alphabet=DNA --input-file=input.csv --tau=0.35 --output-path=results --verbose

This creates a directory called :code:`results` in which all the result files
are written.

.. note::

   All command-line arguments and configuration parameters are documented
   :doc:`here <cmd_args>`.

Tutorial II: Step-by-step processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of using a single call to perform all the processing steps, one can
split these steps into multiple function calls. Here's a code sample:

.. code-block:: python
   :linenos:

   from genonets.cmdl_handler import CmdParser
   from genonets.interface import Genonets

   # Parse the command line argument and create the Genonets object
   gn = Genonets(CmdParser().get_args(), process=False)

   # Use 'gn' to create genotype networks for all genotype sets.
   gn.create()

   # Perform all available analyses on all genotype networks.
   gn.analyze()

   # Write all genotype networks to files in GML format. For a genotype network
   # with two or more components, two files are generated: One corresponds to the
   # entire network with all components, and the other corresponds to the dominant
   # component only.
   gn.save()

   # Save all genotype network level measures to 'Genotype_set_measures.txt'.
   gn.save_network_results()

   # Save all genotype level measures to '<genotypeSetName>_genotype_measures.txt'
   # files. One file per genotype set is generated.
   gn.save_genotype_results()

Tutorial III: Selective processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is also possible to process only a selection of genotype sets from the
input data. Also, it is possible to perform only a selection of available
analysis types. Here's a code sample:

.. code-block:: python
   :linenos:

   from genonets.cmdl_handler import CmdParser
   from genonets.interface import Genonets

   # Parse the command line argument and create the Genonets object
   gn = Genonets(CmdParser().get_args(), process=False)

   # Use 'gn' to create genotype networks for all genotype sets.
   gn.create()

   # Perform only 'Robustness' and 'Evolvability' analyses on just two of
   # the genotype sets available in the input file, i.e., 'Foxa2' and 'Bbx'.
   gn.analyze(["Foxa2", "Bbx"], analyses=[ac.ROBUSTNESS, ac.EVOLVABILITY])

   # Write the given genotype networks to files in GML format.
   # For a genotype network with two or more components, two files are generated:
   # One corresponds to the entire network with all components, and the other
   # corresponds to the dominant component only.
   gn.save(["Foxa2", "Bbx"])

   # Save genotype network level measures for the given genotype sets to
   # 'Genotype_set_measures.txt'.
   gn.save_network_results(["Foxa2", "Bbx"])

   # Save all genotype level measures for the given genotype sets to
   # 'Foxa2_genotype_measures.txt' and 'Bbx_genotype_measures.txt' files.
   gn.save_genotype_results(["Foxa2", "Bbx"])

See :ref:`constants_analysis` and :ref:`api_interface` for more information on
a list of analysis constants and API methods.

.. toctree::
   :hidden:

   input_format
   cmd_args
   concepts_and_terminology
   custom_genetic_code
   result_attributes
   api
