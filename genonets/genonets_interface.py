"""
    genonets_interface
    ~~~~~~~~~~~~~~~~~~

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
    ALL = 0

    # Constructor
    def __init__(self, args):
        # Handle program arguments
        self.cmdArgs = CmdArgs(args)

        # Read file and load input data into memory
        self.inDataDict, self.deltaDict, self.seqToRepDict, self.seqLength = \
            self.buildDataDicts(self.cmdArgs.inFilePath)

        # Get the bit-sequence manipulator object corresponding to the
        # given molecule type.
        self.bitManip = self.getBitManip()

        # Get the NetworkUtils object
        self.netBuilder = NetworkBuilder(self.bitManip)

        # Dictionary: Key=Repertoire, Value=Network. Created when required.
        self.repToNetDict = {}

        # Dictionary: Key=Repertoire, Value=Giant. Created when required.
        # If there is only one component, giant=network
        self.repToGiantDict = {}

        # Reference to the analyzer object
        self.analyzer = None

        # Set the VERBOSE flag
        self.VERBOSE = True if self.cmdArgs.verbose else False

    # ----------------------------------------------------------------------
    #   Public interface
    # ----------------------------------------------------------------------

    # Description: 	Creates genotype networks for the given list of genotype set names.
    # Arguments:
    #   'repertoires': 	List of genotype set names. This argument is optional, and
    #                   if absent, results in the creation of all available
    #                   genotype networks.
    #   'parallel':		Flag to indicate whether or not parallel processing should
    #                   be used. This argument is optional, and is set to 'False'
    #                   by default.
    # Return:		No return value.
    def create(self, repertoires=Gc.ALL, parallel=False):
        if self.VERBOSE:
            print("Creating genotype networks:")

        # If a single string is received, convert it into an iterable
        repertoires = [repertoires] if type(repertoires) == str else repertoires

        # If all repertoires should be considered,
        if repertoires == Gc.ALL:
            # Get a list of all repertoires
            repertoires = self.getRepertoires()

        # If multiprocessing should be used,
        if parallel:
            self.createNets_parallel(repertoires)
        else:
            self.createNets(repertoires)

        if self.VERBOSE:
            sys.stdout.write("Done.\n")

    # Description:	Performs the given set of analyses on the given list of
    #               genotype sets.
    # Arguments:
    #   'repertoires': 	List of genotype set names. This argument is optional, and
    #                   if absent, results in the processing of all available
    #                   genotype networks.
    #   'analyses':		List of analysis types. This argument is optional. If the
    #                   the argument is not provided, all available analyses are
    #                   performed.
    #   'parallel':		Flag to indicate whether or not parallel processing should
    #                   be used. This argument is optional, and is set to 'False'
    #                   by default.
    # Return:		No return value.
    def analyze(self, repertoires=Gc.ALL, analyses=Gc.ALL, parallel=False):
        if self.VERBOSE:
            sys.stdout.write("\nPerforming analyses:")

        # If all repertoires should be considered,
        if repertoires == Gc.ALL:
            # Get a list of all repertoires
            repertoires = self.getRepertoires()

        # If a single string is received, convert it into an iterable
        repertoires = [repertoires] if type(repertoires) == str else repertoires

        # If overlap in one of the requested analyses, there need to be at
        # at least two repertoires in the dataset
        if analyses == Gc.ALL or Ac.OVERLAP in analyses:
            if len(repertoires) < 2:
                print("Error: " +
                      ErrorCodes.getErrDescription(ErrorCodes.NOT_ENOUGH_REPS_OLAP) +
                      ": Tau=" + str(self.cmdArgs.tau))

                raise GenonetsError(
                    ErrorCodes.NOT_ENOUGH_REPS_OLAP,
                    "Tau=" + str(self.cmdArgs.tau))

        # If multiprocessing should be used,
        if parallel:
            # Perform all analyses in parallel; overlap will be ignored.
            self.analyzeNets_parallel(repertoires, analyses)

            if analyses == Gc.ALL or Ac.OVERLAP in analyses:
                # Reset analysis handler to make sure it references
                # the updated dicts
                del self.analyzer
                self.analyzer = AnalysisHandler(self)

                # Use serial processing to perform overlap analysis
                self.analyzeNets(repertoires, [Ac.OVERLAP])
        else:
            # Perform all analyses using serial processing
            self.analyzeNets(repertoires, analyses)

    # Description:	Creates the phenotype network from the given list of
    #               genotype sets.
    # Arguments:
    #   'collection':	Name to give to the phenotype network.
    #   'repertoires': 	List of genotype set names. This argument is optional, and
    #                   if absent, results in the processing of all available
    #                   genotype networks.
    # Return:		igraph object corresponding to the phenotype network.
    # TODO: Add a check to make sure the evolvability analysis has been done
    #       for all networks that are to be processed here ...
    def getEvoNet(self, collection="phenotype_network", repertoires=Gc.ALL):
        # If a single string is received, convert it into an iterable
        repertoires = [repertoires] if type(repertoires) == str else repertoires

        # If all repertoires should be considered,
        if repertoires == Gc.ALL:
            # Get a list of all repertoires
            repertoires = self.getRepertoires()

        # Create a list of giants for the given repertoires
        giants = [self.repToGiantDict[repertoire] for repertoire in repertoires]

        # Create the evolvability network, and get the igraph object
        evoNet = self.netBuilder.createEvoNet(collection, giants)

        return evoNet

    # Description:	Returns the list of names of all genotype sets for which
    #               genotype networks have been created.
    # Return:		List of names of genotype sets.
    def getRepertoires(self):
        repertoires = self.inDataDict.keys()
        repertoires = [repertoires] if type(repertoires) == str else repertoires

        return repertoires

    # Description:	Returns the igraph object for the network corresponding to the
    #               given genotype set name.
    # Arguments:
    #   'repertoire':	Name of the genotype network for which the igraph object
    #                   is required.
    # Return:		igraph object for the required network.
    def getNetworkFor(self, repertoire):
        try:
            return self.repToNetDict[repertoire]
        except KeyError:
            return None

    # Description:	Returns the igraph object for the giant component
    #               corresponding to the given genotype set name.
    # Arguments:
    #   'repertoire':	Name of the genotype network for which the igraph object
    #                   is required.
    # Return:		igraph object for the giant component corresponding to the
    #               given genotype set name.
    def getDominantNetFor(self, repertoire):
        try:
            return self.repToGiantDict[repertoire]
        except KeyError:
            return None

    # Description:	Returns the overlap matrix for all the genotype networks.
    # Return:		Overlap matrix as a list of lists.
    def getOverlapMat(self):
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
            # Perform the overlap compution
            self.analyzer.overlap()

            # Return the resulting matrix
            return self.analyzer.overlapMatrix

    # Description:	Save the given genotype networks as GML files.
    # Arguments:
    #   'repertoires':	List of names of genotype networks to be saved to
    #                   file.
    def save(self, repertoires=Gc.ALL):
        if self.VERBOSE:
            sys.stdout.write("\nWriting GML files for genotype networks ... ")

        # If a single string is received, convert it into an iterable
        repertoires = [repertoires] if type(repertoires) == str else repertoires

        Writer.writeNetsToFile(self.repToNetDict, self.repToGiantDict,
                               self.netBuilder, self.cmdArgs.outPath,
                               WriterFilter.gmlAttribsToIgnore, repertoires)

        if self.VERBOSE:
            sys.stdout.write("Done.\n")

    # Description:	Writes the network level results to file.
    # Arguments:
    #   'repertoires':	List of names of genotype networks for which results
    #                   need to be written to file.
    def saveNetResults(self, repertoires=Gc.ALL):
        if self.VERBOSE:
            sys.stdout.write("\nWriting genotype set level results ... ")

        # If a single string is received, convert it into an iterable
        repertoires = [repertoires] if type(repertoires) == str else repertoires

        Writer.writeNetAttribs(self.repToNetDict, self.repToGiantDict,
                               self.netBuilder, self.cmdArgs.outPath,
                               WriterFilter.netAttribsToIgnore,
                               WriterFilter.net_attribute_to_order,
                               repertoires)

        if self.VERBOSE:
            sys.stdout.write("Done.\n")

    # Description:	Writes the genotype level results to file.
    # Arguments:
    #   'repertoires':	List of names of genotype networks for which results
    #                   need to be written to file.
    def saveGenotypeResults(self, repertoires=Gc.ALL):
        if self.VERBOSE:
            sys.stdout.write("\nWriting genotype level results ... ")

        # If a single string is received, convert it into an iterable
        repertoires = [repertoires] if type(repertoires) == str else repertoires

        Writer.writeSeqAttribs(self.repToGiantDict,
                               self.cmdArgs.outPath,
                               WriterFilter.seqAttribsToIgnore,
                               WriterFilter.seq_attribute_to_order,
                               repertoires)

        if self.VERBOSE:
            sys.stdout.write("Done.\n")

    # Description:	Save the phenotype network as a GML file.
    # Arguments:
    #   'phenotypeNetwork':	igraph object corresponding to the phenotype
    #                       network.
    def saveEvoNet(self, phenotypeNetwork):
        if self.VERBOSE:
            sys.stdout.write("\nWriting GML file for phenotype network ... ")

        Writer.writeNetToFile(phenotypeNetwork, self.cmdArgs.outPath,
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

    def buildDataDicts(self, inFilePath):
        return InReader.buildDataDicts(inFilePath, self.cmdArgs.tau,
                                       self.cmdArgs.moleculeType)

    def getBitManip(self):
        return BitManipFactory.getBitSeqManip(self.cmdArgs.moleculeType,
                                              self.seqLength,
                                              self.cmdArgs.useIndels,
                                              self.cmdArgs.use_reverse_complements)

    def getBitSeqsAndScores(self, repertoire):
        # Get the list of sequences for the given repertoire
        sequences = self.inDataDict[repertoire].keys()

        # Get the list of scores for the given repertoire
        scores = [self.inDataDict[repertoire][seq] for seq in sequences]

        return sequences, scores

    # Create genotype networks for the given, or all repertoires
    def createNets(self, repertoires):
        # For each repertoire,
        for repertoire in repertoires:
            if self.VERBOSE:
                sys.stdout.write(repertoire + " ... ")

            # Get the sequences and scores
            seqs, scores = self.getBitSeqsAndScores(repertoire)

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
    def createNets_parallel(self, repertoires):
        # Instantiate a concurrent queue for results
        resultsQueue = Queue()

        # Compute indices to be used in the loop
        indices = Genonets.process_blocks(len(repertoires), self.cmdArgs.num_procs)

        for i in indices:
            # Create separate processes for each repertoire in the current block
            processes = [
                Process(target=Genonets.createGN,
                        args=(self.inDataDict[repertoires[j]],
                              self.cmdArgs,
                              self.seqLength,
                              resultsQueue,
                              repertoires[j])
                        )
                for j in range(i - 1, Genonets.len_finished_reps(
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
            while len(self.repToNetDict) != Genonets.len_finished_reps(
                    i, len(repertoires), self.cmdArgs.num_procs):
                result = resultsQueue.get()

                if self.VERBOSE:
                    sys.stdout.write(result[0] + " ... ")

                self.repToNetDict[result[0]] = result[1][0]
                self.repToGiantDict[result[0]] = result[1][1]

    @staticmethod
    def createGN(seqScrDict, args, seqLength, resultsQueue, repertoire):
        # Get a reference to the bit manipulator
        bitManip = BitManipFactory.getBitSeqManip(args.moleculeType,
                                                  seqLength, args.useIndels)

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

    def analyzeNets(self, repertoires, analyses):
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

    def analyzeNets_parallel(self, repertoires, analyses):
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
        indices = Genonets.process_blocks(len(repertoires), self.cmdArgs.num_procs)

        for i in indices:
            # Create separate processes for each repertoire in the current block
            processes = [
                Process(target=Genonets.analyzeGN,
                        args=(copy.deepcopy(self_copy),
                              analyses,
                              resultsQueue,
                              repertoires[j])
                        )
                for j in range(i - 1, Genonets.len_finished_reps(
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
            while len(self.repToNetDict) != Genonets.len_finished_reps(
                    i, len(repertoires), self.cmdArgs.num_procs):

                result = resultsQueue.get()

                if self.VERBOSE:
                    print("\tAnalysis results received for: " + result[0])

                self.repToNetDict[result[0]] = result[1][0]
                self.repToGiantDict[result[0]] = result[1][1]

    @staticmethod
    def process_blocks(num_repertoires, num_processes):
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
    def len_finished_reps(cur_index, len_repertoires, max_procs):
        return min(len_repertoires, (cur_index - 1) + max_procs)

    def analyzeNets_parallel_old(self, repertoires, analyses):
        # Instantiate a concurrent queue for results
        resultsQueue = Queue()

        # Create separate processes for each repertoire
        processes = [
            Process(target=Genonets.analyzeGN,
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
    def analyzeGN(genonetsCopy, analyses, resultsQueue, repertoire):
        # Initialize the AnalysisHandler object
        analyzer = AnalysisHandler(genonetsCopy, parallel=True)

        # Perform the analyses
        analyzer.analyze(repertoire, analyses)

        # Create a tuple in which to put the results
        resultTuple = (
            repertoire,
            (
                genonetsCopy.getNetworkFor(repertoire),
                genonetsCopy.getDominantNetFor(repertoire)
            )
        )

        # Add results to the shared queue
        resultsQueue.put(resultTuple)

        # Close the queue for this process
        resultsQueue.close()
