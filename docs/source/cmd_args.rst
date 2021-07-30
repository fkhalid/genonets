Command-line arguments and Configuration parameters
===================================================

Parameters can either be provided on the command-line, or specified in a
configuration file. If a configuration file is used, the only mandatory
command-line argument is :code:`--config-file` (described
:ref:`here <config_file_name>`). Otherwise, all :ref:`mandatory_args` must be
provided.

The following command is an example of using command-line arguments:

.. code-block:: python

    python -u -m genonets.sample.minimal --alphabet=DNA --include-indels --input-file=input.csv --tau=0.35 --output-path=results --verbose

.. _mandatory_args:

Mandatory parameters
^^^^^^^^^^^^^^^^^^^^

.. _config_file_name:

1. :code:`--config-file`: The complete path to the :ref:`config_file`.

   When specified, all other arguments are treated as
   optional. Any other arguments provided in addition to this will override the
   values provided in the configuration file. It is sufficient to provide all
   other arguments in the configuration file, significantly simplifying the
   command-line.

2. :code:`--input-file`: The complete path to the input file, relative to the
   current directory. The file must be in the :doc:`input_format`.

3. :code:`--output-path`: Path to the directory where output files should be
   generated. The directory will be created if it does not already exist.

4. :code:`--alphabet`: Each genotype in the input file must be a string of
   letters from the alphabet specified here. Available options are :code:`RNA`,
   :code:`DNA`, :code:`Protein`, and :code:`Binary`.

.. _tau:

5. :code:`tau`: Minimum score threshold: all genotypes in the input file with
   :ref:`input_format_score` values below this threshold will be ignored. Tau
   must be a number.

.. _optional_args:

Optional parameters
^^^^^^^^^^^^^^^^^^^

1. :code:`--include-indels`: If specified, mutations that shift the entire
   genotype sequence by one letter are also considered. Only point mutations
   are considered otherwise.

   Please note that when used on the command-line, this parameter does not
   require a value; merely specifying it is sufficient.

2. :code:`--use-reverse-complements`: When specified, reverse complements are
   considered during creation and analysis of genotype networks for alphabet
   type DNA. This option is not valid for other alphabet types.

   Please note that when used on the command-line, this parameter does not
   require a value; merely specifying it is sufficient.

.. _args_store_epistasis_squares:

3. :code:`--store-epistasis-squares`: If this parameter is specified and
   epistasis analysis is requested, all computed squares are stored on disk.

   Please note that when used on the command-line, this parameter does not
   require a value; merely specifying it is sufficient.

.. _use_all_components:

4. :code:`--use-all-components`: If specified, evolvability analysis considers
   all connected components, not just the dominant network. It is currently not
   possible to consider all connected components for analyses other than
   evolvability.

   Please note that when used on the command-line, this parameter does not
   require a value; merely specifying it is sufficient.

5. :code:`--num-processes`: No. of processes to be used in parallel processing.
   The value must be a positive integer.

6. :code:`--verbose`: Information about processing steps is written to
   :code:`stdout`.

   Please note that when used on the command-line, this parameter does not
   require a value; merely specifying it is sufficient.

7. :code:`--genetic-code-file`: Complete path to the file which contains the
   mapping from the genetic code to each letter in the selected alphabet. See
   :doc:`custom_genetic_code` for information about the file format.

8. :code:`--codon-alphabet`: Alphabet to use for the codons in the genetic code
   file. The possible values are :code:`DNA` and :code:`RNA`.

9. :code:`--include-indels-for-codons`: When specified, indels are considered
   when checking for mutations in codons (for custom genetic code).

10. :code:`--use-rc-for-codons`: When specified, reverse complements are
    considered when checking for mutations in codons. This is only applicable
    when the alphabet for codons is DNA, and a genetic code file is specified.

.. _config_file:

Configuration file
^^^^^^^^^^^^^^^^^^

In order to simplify the command-line, all command-line arguments can be loaded
from a configuration file. E.g.,

.. code-block:: python

    python -u -m genonets.sample.minimal -c config.cfg

The above command will load all command-line arguments from :code:`config.cfg`.

Following is a sample :code:`config.cfg`:

.. code-block::

    [General]
    input-file: input.csv
    output-path: results/minimal
    alphabet: DNA
    tau: 0.35
    include-indels: yes
    use-reverse-complements: no
    store-epistasis-squares: no
    use-all-components: no
    num-processes: 1
    verbose: yes

    [Genetic Code]
    genetic-code-file: non_standard_code.csv
    codon-alphabet: DNA
    include-indels-for-codons: no
    use-rc-for-codons: no

Please note that the :code:`Genetic Code` section is optional; it can either be
removed or commented out. Placing the :code:`#` character at the beginning of a
line comments it out.