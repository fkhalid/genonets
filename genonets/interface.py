"""
    Public interface to Genonets functions.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import sys
import copy
from multiprocessing import Process, Queue

from genonets.config import ConfigParser
from genonets.writer import Writer
from genonets.reader import InReader
from genonets.reader import GeneticCodeReader
from genonets.graph_utils import NetworkBuilder
from genonets.seq_bit_impl import BitManipFactory
from genonets.filters import WriterFilter
from genonets.analysis_handler import AnalysisHandler
from genonets.constants import ErrorCodes
from genonets.constants import GenonetsConstants as Gc
from genonets.constants import AnalysisConstants as Ac


# TODO: Review and refactor the code for parallel network creation and
#   analysis, as the current code appears outdated ...
class Genonets:
    """
    Encapsulates the Genonets public API.

    """

    ALL = 0

    def __init__(self, args, process=False, parallel=False):
        """
        Initiate parsing of the input file, and load the parsed data into a
        `Genonets` object.

        A simple way to create a `Genonets` object is the following::

            gn = Genonets(CmdParser(args).getArgs())

        where, `CmdParser` can be imported as follows::

            from genonets.cmdl_handler import CmdParser

        The `args` variable is a list of command line arguments. Please see
        the `genonets` package documentation for the list and descriptions of
        all available command line arguments.

        :param args: A populated `CmdArgs` object.
        :param process: If 'True', in addition to creating the object,
                    initiates complete processing, i,e., creates genotype
                    networks for all genotype sets in the input data, performs
                    all available analyses on all genotype networks, and
                    generates all result files.
        :param parallel: Flag to indicate whether or not parallel processing
                    should be used. This parameter is only useful with
                    'process=True'.

        """

        # Handle program arguments
        self.cmd_args = ConfigParser(args)

        # Set the VERBOSE flag
        self.VERBOSE = True if self.cmd_args.verbose else False

        # Dict {letter: list of letters that are 1-neighbors}
        # Should be 'None' if the codon-to-letter map should not be used
        self.letter_to_neighbors = self._create_letter_to_neighbor_map()

        # if self.VERBOSE:
        #     print("\nLoading input data ... ")
        print("\nLoading input data ... ")

        # Read file and load input data into memory
        self.inDataDict,\
            self.deltaDict, \
            self.seqToRepDict, \
            self.seqLength, \
            self.ordered_genotype_sets = self._build_data_dicts(
                self.cmd_args.in_file_path, self.letter_to_neighbors
            )

        if self.VERBOSE:
            print("Done.\n")

        # Pass the ordered list of genotype sets to the WriterFilter
        WriterFilter.ORDERED_GENOTYPE_SETS = self.ordered_genotype_sets

        # Get the bit-sequence manipulator object corresponding to the
        # given molecule type.
        self._alphabet = BitManipFactory.get_alphabet(
            self.cmd_args.molecule_type
        )
        self.bitManip = self._bit_manipulator()

        # Get the NetworkBuilder object
        self.netBuilder = NetworkBuilder(
            self.bitManip,
            self.cmd_args.use_reverse_complements,
            self.VERBOSE
        )

        # Dictionary: Key=Repertoire, Value=Network. Created when required.
        self.repToNetDict = {}

        # Dictionary: Key=Repertoire, Value=Giant. Created when required.
        # If there is only one component, giant=network
        self.repToGiantDict = {}

        self.alternate_bit_manip = None
        self.alternate_net_builder = None
        self.rep_to_alternate_net_dict = None
        self.rep_to_alternate_giant_dict = None

        if self.letter_to_neighbors:
            # Bit manipulator and net builder for the alternate network that
            # supports custom alphabet but does not take genetic code into
            # consideration
            self.alternate_bit_manip = BitManipFactory.getBitSeqManip(
                moleculeType=self.cmd_args.molecule_type,
                seqLength=self.seqLength,
                useIndels=self.cmd_args.use_indels,
                letter_to_neighbors=None,
                alphabet=list(self.letter_to_neighbors.keys()),
                useReverseComplements=False
            )

            self.alternate_net_builder = NetworkBuilder(
                self.alternate_bit_manip,
                False,
                self.VERBOSE
            )

            self.rep_to_alternate_net_dict = {}
            self.rep_to_alternate_giant_dict = {}

        # Phenotype network
        self.pheno_net = None

        # Reference to the analyzer object
        self.analyzer = None

        # If the user has requested complete processing with default settings,
        if process:
            # Perform all processing steps
            self._process_all(parallel)

    def create(self, genotype_sets=Gc.ALL, parallel=False):
        """
        Create genotype networks for the given list of genotype set names.

        :param genotype_sets: List of names of the genotype sets for which the
                    genotype networks should be created. If a value is not
                    explicitly specified for this parameter, genotype networks
                    are constructed for all genotype sets available in the
                    parsed data.
        :param parallel: Flag to indicate whether or not parallel processing
                    should be used.

        :return: No return value

        """

        if self.VERBOSE:
            print("Creating genotype networks:")

        # If a single string is received, convert it into an iterable
        if isinstance(genotype_sets, str):
            genotype_sets = [genotype_sets]

        # If all genotype_sets should be considered,
        if genotype_sets == Gc.ALL:
            # Get a list of all genotype_sets
            genotype_sets = self.genotype_sets()

        # If multiprocessing should be used,
        if parallel:
            self._create_networks_parallel(genotype_sets)

            # If a genetic code is in use
            if self.letter_to_neighbors:
                self._create_alternate_networks_parallel(genotype_sets)
        else:
            self._create_networks(genotype_sets)

        if self.VERBOSE:
            sys.stdout.write("Done.\n")

    def analyze(self, genotype_sets=Gc.ALL, analyses=Gc.ALL, parallel=False):
        """
        Performs all analyses provided in the list of analysis types, on the
        given genotype sets.

        This method can only be used if `create` has already been called on
        the same `Genonets` object.

        :param genotype_sets: List of names of the genotype sets for which the
                    genotype networks should be created. If a value is not
                    explicitly specified for this parameter, genotype networks
                    are constructed for all genotype sets available in the
                    parsed data.
        :param analyses: List of analysis type constants. These constants are
                    defined in the class
                    `genonets.genonets_constants.AnalysisConstants`. If the
                    value for this parameter is not explicitly set, all
                    available analyses are performed.
        :param parallel: Flag to indicate whether or not parallel processing
                    should be used.

        :return: No return value.

        """

        if self.VERBOSE:
            sys.stdout.write("\nPerforming analyses:")

        # If all genotype_sets should be considered,
        if genotype_sets == Gc.ALL:
            # Get a list of all genotype_sets
            genotype_sets = self.genotype_sets()

        # If a single string is received, convert it into an iterable
        if isinstance(genotype_sets, str):
            genotype_sets = [genotype_sets]

        # If overlap in one of the requested analyses, there need to be at
        # at least two genotype_sets in the dataset
        if analyses == Gc.ALL or Ac.OVERLAP in analyses:
            if len(genotype_sets) < 2:
                error_description = ErrorCodes.getErrDescription(
                    ErrorCodes.NOT_ENOUGH_REPS_OLAP
                )
                print(f'Error: {error_description} : Tau={self.cmd_args.tau}')

        # If multiprocessing should be used,
        if parallel:
            # Perform all analyses in parallel; overlap will be ignored.
            self._analyze_networks_parallel(genotype_sets, analyses)

            if analyses == Gc.ALL or Ac.OVERLAP in analyses:
                # Reset analysis handler to make sure it references
                # the updated dicts
                del self.analyzer
                self.analyzer = AnalysisHandler(self)

                # Use serial processing to perform overlap analysis
                self._analyze_networks(genotype_sets, [Ac.OVERLAP])
        else:
            # Perform all analyses using serial processing
            self._analyze_networks(genotype_sets, analyses)

    def phenotype_network(self, collection_name="phenotype_network", genotype_sets=Gc.ALL):
        """
        Create the phenotype network from the given list of genotype sets.

        :param collection_name: The name to be assigned to the phenotype
                    network.
        :param genotype_sets: List of names of the genotype sets for which the
                    phenotype network should be created. If a value is not
                    explicitly specified for this parameter, all available
                    genotype sets are considered.

        :return: `igraph.Graph` object representing the phenotype network.

        """

        # If a single string is received, convert it into an iterable
        if isinstance(genotype_sets, str):
            genotype_sets = [genotype_sets]

        # If all genotype_sets should be considered,
        if genotype_sets == Gc.ALL:
            # Get a list of all genotype_sets
            genotype_sets = self.genotype_sets()

        # Create a list of giants for the given genotype_sets
        giants = [
            self.repToGiantDict[repertoire]
            for repertoire in genotype_sets
        ]

        # Create the phenotype network, and get the igraph object
        self.pheno_net = self.netBuilder.createEvoNet(collection_name, giants)

        return self.pheno_net

    def genotype_sets(self):
        """
        Get a list of names of all genotype sets for which genotype networks
        have been created.

        :return: List of names of genotype sets.

        """

        repertoires = [*self.inDataDict]

        repertoires = [repertoires] if type(repertoires) == str else repertoires

        return repertoires

    def genotype_network(self, genotype_set):
        """
        Get the `igraph` object for the network corresponding to the given
        genotype set name.

        The `igraph` object in this case refers to the entire network, i.e.,
        all connected components.

        Note: This method can only be used if the genotype network
        corresponding to the requested genotype set name has already been
        created.

        :param genotype_set: Name of the genotype set for which the genotype
                    network is requested.

        :return: Object of type `igraph.Graph`.

        """

        try:
            return self.repToNetDict[genotype_set]
        except KeyError:
            return None

    def dominant_network(self, genotype_set):
        """
        Get the `igraph` object for the *dominant* network corresponding to
        the given genotype set name.

        The dominant network refers to the giant component in the network.

        Note: This method can only be used if the genotype network
        corresponding to the requested genotype set name has already been
        created.

        :param genotype_set: Name of the genotype set for which the
                    genotype network is requested.

        :return: Object of type `igraph.Graph`.

        """

        try:
            return self.repToGiantDict[genotype_set]
        except KeyError:
            return None

    def save(self, genotype_sets=Gc.ALL):
        """
        Write the genotype networks corresponding to the given genotype sets
        to file.

        The networks are saved in GML format. For networks with more than one
        components, separate files are generated for the entire network and the
        dominant network.

        Note: This method can be used only after `analyze()` has been called
        on the given genotype sets.

        :param genotype_sets: List of names of genotype sets for which the
                    genotype should be written to file. If a value is not
                    explicitly specified for this parameter, result files are
                    written for all genotype sets.

        :return: No return value.

        """

        if self.VERBOSE:
            print('\nWriting GML files for genotype networks ... ', end=' ... ')

        # If a single string is received, convert it into an iterable
        if isinstance(genotype_sets, str):
            genotype_sets = [genotype_sets]

        Writer.writeNetsToFile(
            self.repToNetDict,
            self.repToGiantDict,
            self.netBuilder,
            self.cmd_args.output_path,
            WriterFilter.gmlAttribsToIgnore,
            genotype_sets
        )

        if self.VERBOSE:
            print('Done.\n')

    def save_network_results(self, genotype_sets=Gc.ALL):
        """
        Write the genotype set level results to file.

        A file named 'Genotype_set_measures.txt' is generated in the output
        directory specified at the time of the `Genonets` object creation.

        Note: This method can be used only after `analyze()` has been called
        on the given genotype sets.

        :param genotype_sets: List of names of genotype sets for which to
                    generate the result files. If a value is not explicitly
                    specified for this parameter, result files are written for
                    all genotype sets.

        :return: No return value.

        """

        if self.VERBOSE:
            print('\nWriting genotype set level results ... ', end=' .. ')

        # If a single string is received, convert it into an iterable
        if isinstance(genotype_sets, str):
            genotype_sets = [genotype_sets]

        Writer.writeNetAttribs(
            self.repToNetDict,
            self.repToGiantDict,
            self.netBuilder,
            self.cmd_args.output_path,
            WriterFilter.netAttribsToIgnore,
            WriterFilter.net_attribute_to_order,
            WriterFilter.genotype_set_to_order,
            genotype_sets
        )

        if self.VERBOSE:
            print('Done.\n')

    def save_genotype_results(self, genotype_sets=Gc.ALL):
        """
        Write the genotype level results to files.

        A results file is generated for each genotype set.

        Note: This method can be used only after `analyze()` has been called
        on the given genotype sets.

        :param genotype_sets: List of names of genotype sets for which to
                    generate the result files. If a value is not explicitly
                    specified for this parameter, result files are written for
                    all genotype sets.

        :return: No return value.

        """

        if self.VERBOSE:
            print('\nWriting genotype level results ... ', end=' ... ')

        # If a single string is received, convert it into an iterable
        if isinstance(genotype_sets, str):
            genotype_sets = [genotype_sets]

        Writer.writeSeqAttribs(
            self.repToGiantDict,
            self.cmd_args.output_path,
            WriterFilter.seqAttribsToIgnore,
            WriterFilter.seq_attribute_to_order,
            genotype_sets
        )

        if self.VERBOSE:
            print('Done.\n')

    def save_phenotype_network(self):
        """
        Write the phenotype network to file in GML format.

        Note: This method can only be used after the phenotype network has been
        created.

        :return: No return value.

        """

        if self.VERBOSE:
            print('\nWriting GML file for phenotype network ... ', end=' ... ')

        Writer.writeNetToFile(
            self.pheno_net,
            self.cmd_args.output_path,
            WriterFilter.gmlAttribsToIgnore
        )

        if self.VERBOSE:
            print('Done.\n')

    # Plots the given network.
    def plot(self, network, layout='auto'):
        self.netBuilder.plotNetwork(
            network,
            layout,
            self.cmd_args.output_path
        )

    # ----------------------------------------------------------------------
    #   Private methods, i.e., those that are not supposed to be part of
    #   public interface
    # ----------------------------------------------------------------------

    def _create_letter_to_neighbor_map(self):
        # Dict {letter: list of letters that are 1-neighbors}
        letter_to_neighbors = None

        # If the genetic code should be consulted
        if self.cmd_args.genetic_code_file and self.cmd_args.molecule_type == 'Protein':
            # Load the codon-to-letter map and also the length of the codon
            self.codon_to_letter_map, self.codon_length = \
                GeneticCodeReader.load_codon_to_letter_map(
                    self.cmd_args.genetic_code_file,
                    self.cmd_args.codon_alphabet
                )

            # Bit manipulator for codons
            codon_bitmanip = BitManipFactory.getBitSeqManip(
                moleculeType=self.cmd_args.codon_alphabet,
                seqLength=self.codon_length,
                useIndels=self.cmd_args.include_indels_for_codons,
                letter_to_neighbors=None,
                useReverseComplements=self.cmd_args.use_rc_for_codons,
                alphabet=BitManipFactory.get_alphabet(
                    self.cmd_args.codon_alphabet)
            )

            # List of codons
            codons = set(self.codon_to_letter_map.keys())

            # Initialize the dict {letter: [list of 1-neighbor letters]}
            letter_to_neighbors = {
                letter: []
                for letter in set(self.codon_to_letter_map.values())
            }

            # For each codon,
            for codon in codons:
                # Create a set with all codons but the current
                other_codons = codons - {codon}

                # Letter encoded by the codon being processed
                codon_encoded_letter = self.codon_to_letter_map[codon]

                # For each of the other codons,
                for other in other_codons:
                    # Check if codon and other are 1-neighbors
                    are_neighbors = codon_bitmanip.areNeighbors(
                        codon_bitmanip.seqToBits(codon),
                        codon_bitmanip.seqToBits(other)
                    )

                    if are_neighbors:
                        # Letter encoded by the other codon
                        other_encoded_letter = self.codon_to_letter_map[other]

                        # To the list of neighboring letters of the letter
                        # encoded by the current codon, append the letter
                        # encoded by the other codon.
                        letter_to_neighbors[codon_encoded_letter].append(
                            other_encoded_letter
                        )

                # Make sure the letter encoded by the current codon does not
                # appear in the list of its neighboring letters
                letter_to_neighbors[codon_encoded_letter] = \
                    list(set(letter_to_neighbors[codon_encoded_letter]) - {codon_encoded_letter})

        return letter_to_neighbors

    def _build_data_dicts(self, inFilePath, letter_to_neighbors):
        return InReader.build_data_dicts(
            inFilePath,
            self.cmd_args.tau,
            self.cmd_args.molecule_type,
            letter_to_neighbors
        )

    def _bit_manipulator(self):
        return BitManipFactory.getBitSeqManip(
            moleculeType=self.cmd_args.molecule_type,
            seqLength=self.seqLength,
            useIndels=self.cmd_args.use_indels,
            letter_to_neighbors=self.letter_to_neighbors,
            useReverseComplements=self.cmd_args.use_reverse_complements,
            alphabet=self._alphabet
        )

    def _bitseqs_and_scores(self, repertoire):
        # Get the list of sequences for the given repertoire
        sequences = [*self.inDataDict[repertoire]]

        # Get the list of scores for the given repertoire
        scores = [self.inDataDict[repertoire][seq] for seq in sequences]

        return sequences, scores

    # Carries out all default processing steps with default settings.
    def _process_all(self, parallel):
        self.create(parallel=parallel)
        self.analyze(parallel=parallel)
        self.save()
        self.save_network_results()
        self.save_genotype_results()
        self.phenotype_network()
        self.save_phenotype_network()

    # Create genotype networks for the given, or all repertoires
    def _create_networks(self, repertoires):
        # For each repertoire,
        for repertoire in repertoires:
            if self.VERBOSE:
                sys.stdout.write(repertoire + " ... ")

            # Get the sequences and scores
            seqs, scores = self._bitseqs_and_scores(repertoire)
            # s_c_tuples = list(zip(seqs, scores))
            # s_c_tuples.sort(key=lambda x: x[0])
            # unzipped = list(zip(*s_c_tuples))
            # seqs = list(unzipped[0])
            # scores = list(unzipped[1])

            # Create the genotype network and store it in a
            # dictionary: Key=Repertoire, Value=Network
            self.repToNetDict[repertoire] = \
                self.netBuilder.createGenoNet(repertoire, seqs, scores)

            net = self.repToNetDict[repertoire]
            if self.VERBOSE:
                sys.stdout.write("[Vertices: {}, Edges: {}]".format(
                    net.vcount(), net.ecount()) + " ... ")

            # Get the number of components in the network
            numComponents = len(self.netBuilder.getComponents(
                self.repToNetDict[repertoire]))

            # If there are more than one components,
            if numComponents > 1:
                # Get the giant component
                giant = self.netBuilder.getGiantComponent(
                    self.repToNetDict[repertoire])

                # Set 'name' attribute for giant
                giant["name"] = repertoire + "_dominant"

                # Reference to the giant component
                self.repToGiantDict[repertoire] = giant
            else:
                # The network and giant component are the same
                self.repToGiantDict[repertoire] = self.repToNetDict[repertoire]

    def create_alternate_network(self, repertoire):
        if self.VERBOSE:
            print(f'Creating alternate genotype network:', end='', flush=True)

        seqs, scores = self._bitseqs_and_scores(repertoire)

        self.rep_to_alternate_net_dict[repertoire] = \
            self.alternate_net_builder.createGenoNet(repertoire, seqs, scores)

        alternate_net = self.rep_to_alternate_net_dict[repertoire]
        if self.VERBOSE:
            sys.stdout.write(
                f'Alternate network [Vertices: {alternate_net.vcount()}, '
                f'Edges: {alternate_net.ecount()}] ... ')

        num_components_alt_net = len(
            self.alternate_net_builder.getComponents(
                self.rep_to_alternate_net_dict[repertoire]))

        if num_components_alt_net > 1:
            alternate_giant = self.alternate_net_builder.getGiantComponent(
                self.rep_to_alternate_net_dict[repertoire]
            )
            self.rep_to_alternate_giant_dict[repertoire] = alternate_giant
        else:
            self.rep_to_alternate_giant_dict[repertoire] = \
                self.rep_to_alternate_net_dict[repertoire]

    # Use multiprocessing to create genotype networks
    # for the given, or all repertoires
    def _create_networks_parallel(self, repertoires):
        # Instantiate a concurrent queue for results
        resultsQueue = Queue()

        # Compute indices to be used in the loop
        indices = Genonets._process_blocks(
            len(repertoires), self.cmd_args.num_processes
        )

        for i in indices:
            # Create separate processes for each repertoire in the current block
            processes = [
                Process(target=Genonets._create_gn,
                        args=(self.inDataDict[repertoires[j]],
                              self.cmd_args,
                              self.seqLength,
                              resultsQueue,
                              repertoires[j])
                        )
                for j in range(i - 1, Genonets._len_finished_reps(
                    i, len(repertoires), self.cmd_args.num_processes))
            ]

            # Start the processes
            for p in processes:
                p.start()

            # Spin lock
            # FIXME: The condition in the loop can result in an infinite
            # iteration if one of the processes does not put results
            # in the queue. This condition should be replaced
            # with one that is reliable ...
            while len(self.repToNetDict) != Genonets._len_finished_reps(
                    i, len(repertoires), self.cmd_args.num_processes):
                result = resultsQueue.get()

                if self.VERBOSE:
                    sys.stdout.write(result[0] + " ... ")

                self.repToNetDict[result[0]] = result[1][0]
                self.repToGiantDict[result[0]] = result[1][1]

    @staticmethod
    def _create_gn(seqScrDict, args, seqLength, resultsQueue, repertoire):
        # Get a reference to the bit manipulator
        bitManip = BitManipFactory.getBitSeqManip(
            moleculeType=args.molecule_type,
            seqLength=seqLength,
            useIndels=args.use_indels,
            letter_to_neighbors=None,   # FIXME: make it variable, not fixed ...
            useReverseComplements=args.use_reverse_complements,
            alphabet=BitManipFactory.get_alphabet(args.molecule_type)
        )

        # Get the sequences for the given repertoire
        sequences = [*seqScrDict]

        # Get the list of scores for the given repertoire
        scores = [seqScrDict[sequence] for sequence in sequences]

        # s_c_tuples = list(zip(sequences, scores))
        # s_c_tuples.sort(key=lambda x: x[0])
        # unzipped = list(zip(*s_c_tuples))
        # sequences = list(unzipped[0])
        # scores = list(unzipped[1])

        # Create a network builder object
        netBuilder = NetworkBuilder(
            bitManip, args.use_reverse_complements, verbose=False)

        # Create the genotype network
        network = netBuilder.createGenoNet(repertoire, sequences, scores)

        # Get the number of components in the network
        numComponents = len(netBuilder.getComponents(network))

        # If there are more than one components,
        if numComponents > 1:
            # Get the giant component
            giant = netBuilder.getGiantComponent(network)

            # Set 'name' attribute for giant
            giant["name"] = repertoire + "_dominant"
        else:
            # The network and giant component are the same
            giant = network

        # Create a tuple in which to put the results
        netTuple = (repertoire, (network, giant))

        # Add the created network objects to the shared queue
        resultsQueue.put(netTuple)

        # Close the queue for this process
        resultsQueue.close()

    def _create_alternate_networks_parallel(self, repertoires):
        # Instantiate a concurrent queue for results
        results_queue = Queue()

        # Compute indices to be used in the loop
        indices = self._process_blocks(
            len(repertoires), self.cmd_args.num_processes
        )

        for i in indices:
            # Create separate processes for each repertoire in the current block
            processes = [
                Process(target=Genonets._create_alternate_gn,
                        args=(self.inDataDict[repertoires[j]],
                              self.cmd_args,
                              self.seqLength,
                              list(self.letter_to_neighbors.keys()),
                              results_queue,
                              repertoires[j])
                        )
                for j in range(i - 1, Genonets._len_finished_reps(
                    i, len(repertoires), self.cmd_args.num_processes))
            ]

            # Start the processes
            for p in processes:
                p.start()

            # Spin lock
            # FIXME: The condition in the loop can result in an infinite
            # iteration if one of the processes does not put results
            # in the queue. This condition should be replaced
            # with one that is reliable ...
            while len(self.rep_to_alternate_net_dict) != Genonets._len_finished_reps(
                    i, len(repertoires), self.cmd_args.num_processes):
                result = results_queue.get()

                if self.VERBOSE:
                    sys.stdout.write(result[0] + " ... ")

                self.rep_to_alternate_net_dict[result[0]] = result[1][0]
                self.rep_to_alternate_giant_dict[result[0]] = result[1][1]

    @staticmethod
    def _create_alternate_gn(sequence_to_score_map, args, sequence_length, alphabet, results_queue, repertoire):
        bit_manip = BitManipFactory.getBitSeqManip(
            moleculeType=args.molecule_type,
            seqLength=sequence_length,
            useIndels=args.use_indels,
            letter_to_neighbors=None,
            alphabet=alphabet,
            useReverseComplements=False
        )

        # Get the sequences for the given repertoire
        sequences = [*sequence_to_score_map]

        # Get the list of scores for the given repertoire
        scores = [sequence_to_score_map[sequence] for sequence in sequences]

        # Create a network builder object
        net_builder = NetworkBuilder(
            bit_manip, use_reverse_complements=False, verbose=False)

        # Create the genotype network
        network = net_builder.createGenoNet(repertoire, sequences, scores)

        # Get the number of components in the network
        num_components = len(net_builder.getComponents(network))

        # If there are more than one components,
        if num_components > 1:
            # Get the giant component
            giant = net_builder.getGiantComponent(network)
        else:
            # The network and giant component are the same
            giant = network

        # Create a tuple in which to put the results
        net_tuple = (repertoire, (network, giant))

        # Add the created network objects to the shared queue
        results_queue.put(net_tuple)

        # Close the queue for this process
        results_queue.close()

    def _analyze_networks(self, repertoires, analyses):
        # Instantiate the analyzer. Initializing here makes it possible
        # to perform actions based on the list of requested analyses.
        self.analyzer = AnalysisHandler(self, analyses)

        # For each repertoire,
        for repertoire in repertoires:
            if self.VERBOSE:
                sys.stdout.write("\n" + repertoire + ":")
                sys.stdout.write("\n\t")

            # Perform the analysis
            self.analyzer.analyze(repertoire, analyses)

        if self.VERBOSE:
            print('')

    def _analyze_networks_parallel(self, repertoires, analyses):
        if self.VERBOSE:
            print('')

        # Make a copy of self
        self_copy = copy.deepcopy(self)

        # Delete the existing net dicts
        del self.repToNetDict
        del self.repToGiantDict

        # Re-initialize the deleted dicts
        self.repToNetDict = {}
        self.repToGiantDict = {}

        # Manage alternate dicts as well if genetic code is in use
        if self.letter_to_neighbors:
            del self.rep_to_alternate_net_dict
            del self.rep_to_alternate_giant_dict

            self.rep_to_alternate_net_dict = {}
            self.rep_to_alternate_giant_dict = {}

        # Instantiate a concurrent queue for results
        resultsQueue = Queue()

        # Compute indices to be used in the loop
        indices = Genonets._process_blocks(
            len(repertoires), self.cmd_args.num_processes
        )

        for i in indices:
            # Create separate processes for each repertoire in the current block
            processes = [
                Process(target=Genonets._analyze_gn,
                        args=(copy.deepcopy(self_copy),
                              analyses,
                              resultsQueue,
                              repertoires[j])
                        )
                for j in range(i - 1, Genonets._len_finished_reps(
                    i, len(repertoires), self.cmd_args.num_processes))
            ]

            # Start the processes
            for p in processes:
                p.start()

            # Spin lock
            # FIXME: The condition in the loop can result in an infinite
            # iteration if one of the processes does not put results
            # in the queue. This condition should be replaced
            # with one that is reliable ...
            while len(self.repToNetDict) != Genonets._len_finished_reps(
                    i, len(repertoires), self.cmd_args.num_processes):

                result = resultsQueue.get()

                # if self.VERBOSE:
                #     print("\tAnalysis results received for: " + result[0])
                print("\tAnalysis results received for: " + result[0])

                self.repToNetDict[result[0]] = result[1][0]
                self.repToGiantDict[result[0]] = result[1][1]

                if self.letter_to_neighbors:
                    self.rep_to_alternate_net_dict[result[0]] = result[1][2]
                    self.rep_to_alternate_giant_dict[result[0]] = result[1][3]

    @staticmethod
    def _process_blocks(num_repertoires, num_processes):
        # If there are enough processes to process all repertoires
        # in parallel,
        if num_repertoires <= num_processes:
            # Only one iteration is required
            indices = [1]
        else:
            # Calculate the number of iteration that consume all
            # processes
            num_full_iters = num_repertoires // num_processes

            # Calculate the increment required in the last iteration
            incr_last = num_repertoires % num_processes

            # Create a list of increment values per iteration
            index_incrs = [num_processes for _ in range(num_full_iters)] \
                if num_full_iters != 0 else []

            # If required, add the last increment
            if incr_last != 0:
                index_incrs.append(incr_last)

            # Initialize the list of indices
            indices = [1]

            # Indices from index 1 onwards are calculated using prefix sum
            # of the increment list
            for i in range(0, len(index_incrs) - 1):
                indices.append(index_incrs[i] + indices[i])

        return indices

    @staticmethod
    def _len_finished_reps(cur_index, len_repertoires, max_procs):
        return min(len_repertoires, (cur_index - 1) + max_procs)

    def analyzeNets_parallel_old(self, repertoires, analyses):
        # Instantiate a concurrent queue for results
        resultsQueue = Queue()

        # Create separate processes for each repertoire
        processes = [
            Process(target=Genonets._analyze_gn,
                    args=(copy.deepcopy(self),
                          analyses,
                          resultsQueue,
                          repertoire)
                    )
            for repertoire in repertoires
        ]

        # Start the processes
        for p in processes:
            p.start()

        # Delete the existing net dicts
        del self.repToNetDict
        del self.repToGiantDict

        # Re-initialize the deleted dicts
        self.repToNetDict = {}
        self.repToGiantDict = {}

        # Spin lock
        # FIXME: The condition in the loop can result in an infinite
        #        iteration if one of the processes does not put results
        #        in the queue. This condition should be replaced
        #        with one that is reliable ...
        while len(self.repToNetDict) != len(repertoires):
            result = resultsQueue.get()

            print("Analysis results received for: " + result[0])

            self.repToNetDict[result[0]] = result[1][0]
            self.repToGiantDict[result[0]] = result[1][1]

    @staticmethod
    def _analyze_gn(genonets_copy, analyses, results_queue, repertoire):
        # Initialize the AnalysisHandler object
        analyzer = AnalysisHandler(genonets_copy, parallel=True)

        # Perform the analyses
        analyzer.analyze(repertoire, analyses)

        # Create a tuple in which to put the results
        result_tuple = (
            repertoire,
            (
                genonets_copy.genotype_network(repertoire),
                genonets_copy.dominant_network(repertoire),
                genonets_copy.rep_to_alternate_net_dict,
                genonets_copy.rep_to_alternate_giant_dict
            )
        )

        # Add results to the shared queue
        results_queue.put(result_tuple)

        # Close the queue for this process
        results_queue.close()

    # Description:	Returns the overlap matrix for all the genotype networks.
    # Return:		Overlap matrix as a list of lists.
    def _overlap_matrix(self):
        # Overlap matrix can only be computed if the networks have already
        # been created.
        if len(self.repToGiantDict) == 0:
            # Networks have not been created. Therefore, overlap matrix
            # cannot be computed.
            print("Overlap matrix cannot be computed before network creation.")

            return None

        # If the overlap matrix has already been computed,
        if self.analyzer.overlapMatrix:
            # No need to compute again, just return the existing matrix
            return self.analyzer.overlapMatrix
        else:
            # Perform the overlap computation
            self.analyzer.overlap()

            # Return the resulting matrix
            return self.analyzer.overlapMatrix
