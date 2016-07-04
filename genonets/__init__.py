#!/usr/bin/env python

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

    Tutorial I: The one-liner
    -------------------------

    The following code snippet shows the simplest possible way of using Genonets::

        from genonets.cmdl_handler import CmdParser  # For parsing command line arguments
        from genonets.genonets_interface import Genonets  # Interface to Genonets API

        Genonets(CmdParser().getArgs(), process=True)

    Yes, that's it. Two import statements, and just one line of code to create, analyze, and save, all the
    genotype networks in the input file, as well the phenotype network. In fact, these are the contents of the
    `genonets_exmpl_minimal.py` sample file included in the package.

    Assuming the you have downloaded the sample code and changed directory to 'sample/', we can run
    `genonets_exmpl_minimal.py` from the command line as follows::

        python genonets_exmpl_minimal.py DNA true data/genonets_sample_input.txt 0.35 results_simple

    The command line arguments specified above are all mandatory positional arguments, i.e., one must specify each
    one of these arguments, and in the correct order. Here's the ordered list of arguments and the corresponding
    descriptions:

    1. **Alphabet type:** The type of alphabet used in the input file. Supported values are:
     * RNA
     * DNA
     * Protein
     * Binary
    2. **Include indels:** Whether or not indels should be considered mutation types.
    3. **Input file name:** Path to, and name of the input file.
    4. **Tau:** The minimum score value to consider when reading genotypes from the input file.
    . **Result directory:** Path to the directory in which the result files should be created.

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


    Tutorial II: Step by step processing
    ------------------------------------

    Instead of using a single call to perform all the processing steps, one can split these steps into multiple
    function calls. Here's a code sample::

        # Parse the command line argument
        args = CmdParser().getArgs()

        # Create the Genonets object. This will load the input file into memory.
        gn = Genonets(args)

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

    Tutorial II: Parallel processing
    --------------------------------

    Parallel processing an be used independently in network creation and network analysis. Here's how parallel
    processing can be enabled::

        # Use 'gn' to create genotype networks for all genotype sets in parallel.
        gn.create(parallel=True)

        # Perform all available analyses on all genotype networks in parallel.
        gn.analyze(parallel= True)

    Tutorial III: Selective processing
    ----------------------------------

    It is also possible to process only a selection of genotype sets from the input data. Also, it is possible to
    perform only a selection of available analysis types. Here's code sample::

        # Parse the command line argument
        args = CmdParser().getArgs()

        # Create the Genonets object. This will load the input file into
        # memory.
        gn = Genonets(args)

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

    Tutorial III: Customizing results
    ---------------------------------

    This tutorial illustrates the process of customizing the output by adding information to the result file that would
    not be added by Genonets by default.

    The 'Peaks' analysis is used as an example. By default, the 'Peaks' analysis stores results in a dictionary of the
    form::

        {key=peakID : value=[genotypes in the peak]}.

    The sample code that follows customizes this dictionary by adding the score value corresponding to each
    genotype in the list. The resulting dictionary is of the format::

        {key=peakID : value=[(genotype1, score1), ..., (genotypeN, scoreN)]}

    i.e., it is a list of tuples.::

        # Parse the command line argument
        args = CmdParser().getArgs()

        # Create the Genonets object. This will load the input file into
        # memory.
        gn = Genonets(args)

        # Use 'gn' to create genotype networks for all genotype sets.
        gn.create()

        # Perform 'Peaks' analysis
        gn.analyze(analyses=[ac.PEAKS])

        # At this point, the analysis is done. We now need to extract the
        # data we need, i.e., the peaks dictionary.

        # For each genotype set,
        for genotypeSet in gn.genotype_sets():
            # Get the igraph object for the giant
            giant = gn.dominant_network(genotypeSet)

            # Get the dict of peaks {key=peakId : value=[list of sequences in the peak]}
            peaks = giant["Peaks"]

            # Update the dict of peaks, so that instead of just a list of
            # genotypes, we have a list of tuples (sequence, escore).
            newPeaks = {}

            # For each peak,
            for peak in peaks:
                # Initialize the list of tuples
                seqScrTuples = []

                # For each sequence in this peak,
                for sequence in peaks[peak]:
                    # Find the corresponding vertex in the giant
                    try:
                        vertex = giant.vs.find(sequences=sequence)
                    except ValueError:
                        print("Oops! can't find " + sequence + " in giant.")

                    # Get the escore
                    score = giant.vs[vertex.index]["escores"]

                    # Add the tuple to the list of tuples
                    seqScrTuples.append((sequence, score))

                # Add the peak and the corresponding list of tuples to the
                # new peaks dict
                newPeaks[peak] = seqScrTuples

            # Replace the peaks dict in giant with the new dict. This is
            # useful because now if you use the genonets functions below to
            # save the network and results, the updated peaks dict with tuples
            # will be written automatically to file.
            giant["Peaks"] = newPeaks

        # Save networks to file in GML format
        gn.save()

        # Save the results to file from network level analysis
        gn.save_network_results()


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
