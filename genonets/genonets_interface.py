"""
    Public interface to Genonets functions.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import sys
import copy
from multiprocessing import Process, Queue

from cmdl_handler import CmdArgs
from genonets_writer import Writer
from genonets_reader import InReader
from graph_utils import NetworkBuilder
from seq_bit_impl import BitManipFactory
from genonets_filters import WriterFilter
from analysis_handler import AnalysisHandler
from genonets_constants import ErrorCodes
from genonets_exceptions import GenonetsError
from genonets_constants import GenonetsConstants as Gc
from genonets_constants import AnalysisConstants as Ac


class Genonets:
    """
    Encapsulates the Genonets public API.
    """

    ALL = 0

    def __init__(self, args, process=False, parallel=False):
        """
        Initiate parsing of the input file, and load the parsed data into a `Genonets` object.

        A simple way to create a `Genonets` object is the following::

            gn = Genonets(CmdParser(args).getArgs())

        where, `CmdParser` can be imported as follows::

            from genonets.cmdl_handler import CmdParser

        The `args` variable is a list of command line arguments. Please see the `genonets` package
        documentation for the list and descriptions of all available command line arguments.

        :param args: A populated `CmdArgs` object.
        :param process: If 'True', in addition to creating the object, initiates complete processing, i,e., creates
                        genotype networks for all genotype sets in the input data, performs all available analyses on
                        all genotype networks, and generates all result files.
        :param parallel: Flag to indicate whether or not parallel processing should be used. This parameter is only
                         useful with 'process=True'.
        """

        # Handle program arguments
        self.cmdArgs = CmdArgs(args)

        # Read file and load input data into memory
        self.inDataDict, self.deltaDict, self.seqToRepDict, self.seqLength = \
            self._build_data_dicts(self.cmdArgs.inFilePath)

        # Get the bit-sequence manipulator object corresponding to the
        # given molecule type.
        self.bitManip = self._bit_manipulator()

        # Get the NetworkUtils object
        self.netBuilder = NetworkBuilder(self.bitManip)

        # Dictionary: Key=Repertoire, Value=Network. Created when required.
        self.repToNetDict = {}

        # Dictionary: Key=Repertoire, Value=Giant. Created when required.
        # If there is only one component, giant=network
        self.repToGiantDict = {}

        # Phenotype network
        self.pheno_net = None

        # Reference to the analyzer object
        self.analyzer = None

        # Set the VERBOSE flag
        self.VERBOSE = True if self.cmdArgs.verbose else False

        # If the user has requested complete processing with default settings,
        if process:
            # Perform all processing steps
            self._process_all(parallel)

    def create(self, genotype_sets=Gc.ALL, parallel=False):
        """
        Create genotype networks for the given list of genotype set names.

        :param genotype_sets: List of names of the genotype sets for which the genotype
                             networks should be created. If a value is not explicitly
                             specified for this parameter, genotype networks are
                             constructed for all genotype sets available in the parsed
                             data.
        :param parallel: Flag to indicate whether or not parallel processing should
                         be used.
        :return: No return value
        """

        if self.VERBOSE:
            print("Creating genotype networks:")

        # If a single string is received, convert it into an iterable
        genotype_sets = [genotype_sets] if type(genotype_sets) == str else genotype_sets

        # If all genotype_sets should be considered,
        if genotype_sets == Gc.ALL:
            # Get a list of all genotype_sets
            genotype_sets = self.genotype_sets()

        # If multiprocessing should be used,
        if parallel:
            self._create_networks_parallel(genotype_sets)
        else:
            self._create_networks(genotype_sets)

        if self.VERBOSE:
            sys.stdout.write("Done.\n")

    def analyze(self, genotype_sets=Gc.ALL, analyses=Gc.ALL, parallel=False):
        """
        Performs all analyses provided in the list of analysis types, on the given genotype sets.

        This method can only be used if `create` has already been called on the same `Genonets`
        object.

        :param genotype_sets: List of names of the genotype sets for which the genotype
                            networks should be created. If a value is not explicitly
                            specified for this parameter, genotype networks are
                            constructed for all genotype sets available in the parsed
                            data.
        :param analyses: List of analysis type constants. These constants are defined in the class
                         `genonets.genonets_constants.AnalysisConstants`. If the value for this
                         parameter is not explicitly set, all available analyses are performed.
        :param parallel: Flag to indicate whether or not parallel processing should
                         be used.
        :return: No return value.
        """

        if self.VERBOSE:
            sys.stdout.write("\nPerforming analyses:")

        # If all genotype_sets should be considered,
        if genotype_sets == Gc.ALL:
            # Get a list of all genotype_sets
            genotype_sets = self.genotype_sets()

        # If a single string is received, convert it into an iterable
        genotype_sets = [genotype_sets] if type(genotype_sets) == str else genotype_sets

        # If overlap in one of the requested analyses, there need to be at
        # at least two genotype_sets in the dataset
        if analyses == Gc.ALL or Ac.OVERLAP in analyses:
            if len(genotype_sets) < 2:
                print("Error: " +
                      ErrorCodes.getErrDescription(ErrorCodes.NOT_ENOUGH_REPS_OLAP) +
                      ": Tau=" + str(self.cmdArgs.tau))

                raise GenonetsError(
                    ErrorCodes.NOT_ENOUGH_REPS_OLAP,
                    "Tau=" + str(self.cmdArgs.tau))

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

        :param collection_name: The name to be assigned to the phenotype network.
        :param genotype_sets: List of names of the genotype sets for which the phenotype
                              network should be created. If a value is not explicitly
                              specified for this parameter, all available genotype sets
                              are considered.
        :return: `igraph.Graph` object representing the phenotype network.
        """

        # If a single string is received, convert it into an iterable
        genotype_sets = [genotype_sets] if type(genotype_sets) == str else genotype_sets

        # If all genotype_sets should be considered,
        if genotype_sets == Gc.ALL:
            # Get a list of all genotype_sets
            genotype_sets = self.genotype_sets()

        # Create a list of giants for the given genotype_sets
        giants = [self.repToGiantDict[repertoire] for repertoire in genotype_sets]

        # Create the phenotype network, and get the igraph object
        self.pheno_net = self.netBuilder.createEvoNet(collection_name, giants)

        return self.pheno_net

    def genotype_sets(self):
        """
        Get a list of names of all genotype sets for which genotype networks have been created.

        :return: List of names of genotype sets.
        """
        repertoires = self.inDataDict.keys()
        repertoires = [repertoires] if type(repertoires) == str else repertoires

        return repertoires

    def genotype_network(self, genotype_set):
        """target="_blank"
        Get the `igraph` object for the network corresponding to the given genotype set name.

        The `igraph` object in this case refers to the entire network, i.e., all connected
        components.

        Note: This method can only be used if the genotype network corresponding to the requested
        genotype set name has already been created.

        :param genotype_set: Name of the genotype set for which the genotype network is
                             requested.
        :return: Object of type `igraph.Graph`.
        """

        try:
            return self.repToNetDict[genotype_set]
        except KeyError:
            return None

    def dominant_network(self, genotype_set):
        """
        Get the `igraph` object for the *dominant* network corresponding to the given genotype set name.

        The dominant network refers to the giant component in the network.

        Note: This method can only be used if the genotype network corresponding to the requested
        genotype set name has already been created.

        :param genotype_set: Name of the genotype set for which the genotype network is
                             requested.
        :return: Object of type `igraph.Graph`.
        """

        try:
            return self.repToGiantDict[genotype_set]
        except KeyError:
            return None

    def save(self, genotype_sets=Gc.ALL):
        """
        Write the genotype networks corresponding to the given genotype sets to file.

        The networks are saved in GML format. For networks with more than one
        components, separate files are generated for the entire network and the
        dominant network.

        Note: This method can be used only after `analyze()` has been called on the
        given genotype sets.

        :param genotype_sets: List of names of genotype sets for which the genotype
                              should be written to file. If a value is not explicitly specified
                              for this parameter, result files are written for all
                              genotype sets.
        :return: No return value.
        """

        if self.VERBOSE:
            sys.stdout.write("\nWriting GML files for genotype networks ... ")

        # If a single string is received, convert it into an iterable
        genotype_sets = [genotype_sets] if type(genotype_sets) == str else genotype_sets

        Writer.writeNetsToFile(self.repToNetDict, self.repToGiantDict,
                               self.netBuilder, self.cmdArgs.outPath,
                               WriterFilter.gmlAttribsToIgnore, genotype_sets)

        if self.VERBOSE:
            sys.stdout.write("Done.\n")

    def save_network_results(self, genotype_sets=Gc.ALL):
        """
        Write the genotype set level results to file.

        A file named 'Genotype_set_measures.txt' is generated in the output directory
        specified at the time of the `Genonets` object creation.

        Note: This method can be used only after `analyze()` has been called on the
        given genotype sets.

        :param genotype_sets: List of names of genotype sets for which to generate the
                              result files. If a value is not explicitly specified
                              for this parameter, result files are written for all
                              genotype sets.
        :return: No return value.
        """

        if self.VERBOSE:
            sys.stdout.write("\nWriting genotype set level results ... ")

        # If a single string is received, convert it into an iterable
        genotype_sets = [genotype_sets] if type(genotype_sets) == str else genotype_sets

        Writer.writeNetAttribs(self.repToNetDict, self.repToGiantDict,
                               self.netBuilder, self.cmdArgs.outPath,
                               WriterFilter.netAttribsToIgnore,
                               WriterFilter.net_attribute_to_order,
                               genotype_sets)

        if self.VERBOSE:
            sys.stdout.write("Done.\n")

    def save_genotype_results(self, genotype_sets=Gc.ALL):
        """
        Write the genotype level results to files.

        A results file is generated for each genotype set.

        Note: This method can be used only after `analyze()` has been called on the
        given genotype sets.

        :param genotype_sets: List of names of genotype sets for which to generate the
                              result files. If a value is not explicitly specified
                              for this parameter, result files are written for all
                              genotype sets.
        :return: No return value.
        """

        if self.VERBOSE:
            sys.stdout.write("\nWriting genotype level results ... ")

        # If a single string is received, convert it into an iterable
        genotype_sets = [genotype_sets] if type(genotype_sets) == str else genotype_sets

        Writer.writeSeqAttribs(self.repToGiantDict,
                               self.cmdArgs.outPath,
                               WriterFilter.seqAttribsToIgnore,
                               WriterFilter.seq_attribute_to_order,
                               genotype_sets)

        if self.VERBOSE:
            sys.stdout.write("Done.\n")

    def save_phenotype_network(self):
        """
        Write the phenotype network to file in GML format.

        Note: This method can only be used after the phenotype network has been created.

        :return: No return value.
        """

        if self.VERBOSE:
            sys.stdout.write("\nWriting GML file for phenotype network ... ")

        Writer.writeNetToFile(self.pheno_net, self.cmdArgs.outPath,
                              WriterFilter.gmlAttribsToIgnore)

        if self.VERBOSE:
            sys.stdout.write("Done.\n")

    # Plots the given network.
    def plot(self, network, layout="auto"):
        self.netBuilder.plotNetwork(network, layout, self.cmdArgs.outPath)

    # ----------------------------------------------------------------------
    #   Private methods, i.e., those that are not supposed to be part of
    #   public interface
    # ----------------------------------------------------------------------

    def _build_data_dicts(self, inFilePath):
        return InReader.buildDataDicts(inFilePath, self.cmdArgs.tau,
                                       self.cmdArgs.moleculeType)

    def _bit_manipulator(self):
        return BitManipFactory.getBitSeqManip(self.cmdArgs.moleculeType,
                                              self.seqLength,
                                              self.cmdArgs.useIndels,
                                              self.cmdArgs.use_reverse_complements)

    def _bitseqs_and_scores(self, repertoire):
        # Get the list of sequences for the given repertoire
        sequences = self.inDataDict[repertoire].keys()

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

            # Create the genotype network and store it in a
            # dictionary: Key=Repertoire, Value=Network
            self.repToNetDict[repertoire] = \
                self.netBuilder.createGenoNet(repertoire, seqs, scores)

            # Get the number of components in the network
            numComponents = len(self.netBuilder.getComponents(
                self.repToNetDict[repertoire]))

            # If there are more than one components,
            if numComponents > 1:
                # Get the giant component
                giant = self.netBuilder.getGiantComponent(self.repToNetDict[repertoire])

                # Set 'name' attribute for giant
                giant["name"] = repertoire + "_dominant"

                # Reference to the giant component
                self.repToGiantDict[repertoire] = giant
            else:
                # The network and giant component are the same
                self.repToGiantDict[repertoire] = self.repToNetDict[repertoire]

    # Use multiprocessing to create genotype networks
    # for the given, or all repertoires
    def _create_networks_parallel(self, repertoires):
        # Instantiate a concurrent queue for results
        resultsQueue = Queue()

        # Compute indices to be used in the loop
        indices = Genonets._process_blocks(len(repertoires), self.cmdArgs.num_procs)

        for i in indices:
            # Create separate processes for each repertoire in the current block
            processes = [
                Process(target=Genonets._create_gn,
                        args=(self.inDataDict[repertoires[j]],
                              self.cmdArgs,
                              self.seqLength,
                              resultsQueue,
                              repertoires[j])
                        )
                for j in range(i - 1, Genonets._len_finished_reps(
                    i, len(repertoires), self.cmdArgs.num_procs))
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
                    i, len(repertoires), self.cmdArgs.num_procs):
                result = resultsQueue.get()

                if self.VERBOSE:
                    sys.stdout.write(result[0] + " ... ")

                self.repToNetDict[result[0]] = result[1][0]
                self.repToGiantDict[result[0]] = result[1][1]

    @staticmethod
    def _create_gn(seqScrDict, args, seqLength, resultsQueue, repertoire):
        # Get a reference to the bit manipulator
        bitManip = BitManipFactory.getBitSeqManip(args.moleculeType,
                                                  seqLength,
                                                  args.useIndels,
                                                  args.use_reverse_complements)

        # Get the sequences for the given repertoire
        sequences = seqScrDict.keys()

        # Get the list of scores for the given repertoire
        scores = [seqScrDict[sequence] for sequence in sequences]

        # Create a network builder object
        netBuilder = NetworkBuilder(bitManip)

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
            print

    def _analyze_networks_parallel(self, repertoires, analyses):
        if self.VERBOSE:
            print

        # Make a copy of self
        self_copy = copy.deepcopy(self)

        # Delete the existing net dicts
        del self.repToNetDict
        del self.repToGiantDict

        # Re-initialize the deleted dicts
        self.repToNetDict = {}
        self.repToGiantDict = {}

        # Instantiate a concurrent queue for results
        resultsQueue = Queue()

        # Compute indices to be used in the loop
        indices = Genonets._process_blocks(len(repertoires), self.cmdArgs.num_procs)

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
                    i, len(repertoires), self.cmdArgs.num_procs))
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
                    i, len(repertoires), self.cmdArgs.num_procs):

                result = resultsQueue.get()

                if self.VERBOSE:
                    print("\tAnalysis results received for: " + result[0])

                self.repToNetDict[result[0]] = result[1][0]
                self.repToGiantDict[result[0]] = result[1][1]

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
            num_full_iters = num_repertoires / num_processes

            # Calculate the increment required in the last iteration
            incr_last = num_repertoires % num_processes

            # Create a list of increment values per iteration
            index_incrs = [num_processes for i in range(num_full_iters)] \
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
    def _analyze_gn(genonetsCopy, analyses, resultsQueue, repertoire):
        # Initialize the AnalysisHandler object
        analyzer = AnalysisHandler(genonetsCopy, parallel=True)

        # Perform the analyses
        analyzer.analyze(repertoire, analyses)

        # Create a tuple in which to put the results
        resultTuple = (
            repertoire,
            (
                genonetsCopy.genotype_network(repertoire),
                genonetsCopy.dominant_network(repertoire)
            )
        )

        # Add results to the shared queue
        resultsQueue.put(resultTuple)

        # Close the queue for this process
        resultsQueue.close()

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
