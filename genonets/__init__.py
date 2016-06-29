#!/usr/bin/env python

# TODO: Add the download samples directory option, and update code samples accordingly ...

"""
    ========
    Genonets
    ========

    This package provides a high level interface for construction and analysis of genotype networks from data.

    Following features are available:

    * Parsing of genotype-phenotype maps from an input file, provided in the genonets input file format
    * Creation of genotype networks from the input data
    * Various analyses on the constructed genotype networks
    * Generation of result files with attributes from genotype network level analyses, as well as genotype level
      analyses
    * Creation of a phenotype network that shows evolvability and accessibility relationships between the genotype sets
    * Generation of GML files corresponding to the created genotype networks and the phenotype network

    Tutorial
    --------

    The following code snippet shows the simplest possible way of using Genonets::

        from genonets.cmdl_handler import CmdParser  # For parsing command line arguments
        from genonets.genonets_interface import Genonets  # Interface to Genonets API

        Genonets(CmdParser().getArgs(), process=True)

    Yes, that's it. Two import statements, and just one line of code to create, analyze, and save, all the
    genotype networks in the input file, as well the phenotype network. In fact, these are the contents of the
    `genonets_exmpl_minimal.py` sample file included in the package.

    Assuming the `genonets_sample_input.txt` input file is available in the current directory, we can run
    `genonets_exmpl_minimal.py` from the command line as follows::

        python genonets_exmpl_minimal.py DNA true genonets_sample_input.txt 0.35 results_simple

    The command line arguments specified above are all mandatory positional arguments, i.e., one must specify each
    one of these arguments, and in the correct order. Here's the ordered list of arguments and the corresponding
    descriptions:

    #. **Alphabet type:** The type of alphabet used in the input file. Supported values are:
     * RNA
     * DNA
     * Protein
     * Binary
    #. **Include indels:** Whether or not indels should be considered mutation types.
    #. **Input file name:** Path to, and name of the input file.
    #. **Tau:** The minimum score value to consider when reading genotypes from the input file.
    #. **Result directory:** Path to the directory in which the result files should be created.

    Note: These parameters are described in detail in the *Genonets Server Tutorial*.

    In addition to the above listed mandatory arguments, the following optional arguments are also available:

    * **'-rc' or '---use_reverse_complements':** This argument can be specified without a value to enable consideration
      of reverse complements during the construction and analysis of genotype networks. This argument can only be used
      with alphabet type *DNA*.
    * **'-np' or '---num_processes':** The number of processes to use when working in the parallel processing mode. An
      integer value greater than 0 has to be specified as the value to this argument.
    * **'-v' or '---verbose':** This argument can be specified without a value to enable detailed printing of the
      processing steps. Please note that this argument should be used in conjunction with the '-u' option for the
      'python command', e.g.,::

        python -u genonets_exmpl_minimal.py DNA true genonets_sample_input.txt 0.35 results_simple -v



    :author:    Fahad Khalid
    :license:   The MIT License (MIT)

                Copyright (c) 2016 Fahad Khalid and Joshua L. Payne

                Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
                associated documentation files (the "Software"), to deal in the Software without restriction, including
                without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
                copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
                the following conditions:

                The above copyright notice and this permission notice shall be included in all copies or substantial
                portions of the Software.

                THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
                LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
                NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
                WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
                SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
