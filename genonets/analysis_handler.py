"""
    analysis_handler
    ~~~~~~~~~~~~~~~~

    Exposes wrapper functions, one per analysis type. Serves as a collective
    interface to handler classes for all analysis types.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import sys
import json  # For proper list stringification

from genonets_writer import Writer
from path_functions import PathAnalyzer
from landscape_functions import Landscape
from overlap_functions import OverlapAnalyzer
from structure_functions import StructureAnalyzer
from robustness_functions import RobustnessAnalyzer
from evolvability_functions import EvolvabilityAnalyzer
from accessibility_functions import AccessibilityAnalyzer
from genonets_constants import AnalysisConstants as Ac
from genonets_constants import EpistasisConstants as Epi


class AnalysisHandler:
    # Constructor
    def __init__(self, caller, analyses=Ac.ALL, parallel=False):
        # Store a reference to the caller object
        self.caller = caller

        # Set the flag that indicates whether or not this instance
        # is being used for parallel processing
        self.parallel = parallel

        # Dictionary: Key=Repertoire, Value=Network. Created when required.
        self.repToNetDict = self.caller.repToNetDict

        # Dictionary: Key=Repertoire, Value=Giant. Created when required.
        self.repToGiantDict = self.caller.repToGiantDict

        # Reference to the NetBuilder object in use
        self.netBuilder = self.caller.netBuilder

        # Get the bit-sequence manipulator object corresponding to the
        # given molecule type.
        self.bitManip = self.caller.bitManip

        # Reference to input data dictionary
        self.inDataDict = self.caller.inDataDict

        # Reference to delta dict
        self.deltaDict = self.caller.deltaDict

        # Reference to
        self.seqToRepDict = self.caller.seqToRepDict

        # Verbosity flag
        self.VERBOSE = self.caller.VERBOSE

        # Dictionary to store 'analysis type' to 'function'
        # mapping
        self.analysisToFunc = {
            Ac.PEAKS: self.peaks,
            Ac.PATHS: self.paths,
            Ac.PATHS_RATIOS: self.paths_ratios,
            Ac.EPISTASIS: self.epistasis,
            Ac.ROBUSTNESS: self.robustness,
            Ac.EVOLVABILITY: self.evolvability,
            Ac.ACCESSIBILITY: self.accessibility,
            Ac.NEIGHBOR_ABUNDANCE: self.neighborAbundance,
            Ac.PHENOTYPIC_DIVERSITY: self.phenotypicDiversity,
            Ac.STRUCTURE: self.structure,
            Ac.OVERLAP: self.overlap
        }

        # Flag to indicate whether or not the genotypes should be
        # considered double stranded, i.e., whether or not
        # reverse complements should be used in evolvability
        # computations.
        self.isDoubleStranded = self.caller.cmdArgs.use_reverse_complements

        # If 'Evolvability' analysis has been requested, initialize
        # data structures specific to 'Evolvability' analysis
        # Note: This is by design, since building these data
        # structures once is a lot more efficient than building them
        # again and again for each repertoire.
        if analyses == Ac.ALL or Ac.EVOLVABILITY in analyses:
            # Dict - {sequence: [repertoires]}, with only those repertoires
            # for which the sequence in the giant.
            self.seqToRepDict_evo = None

            # Dict - {reverse complement in bit format: original sequence in
            # string format}
            self.rcToSeqDict = None

            # Dict - {sequence in bit format: sequence in string format}
            self.bitsToSeqDict = None

            # Initialize analysis specific data structures used in
            # in the analysis of all repertoires
            self.init_evolvability()

        # Reference to the overlap matrix
        self.overlapMatrix = None

    def getFuncsFor(self, analysis):
        try:
            funcs = self.analysisToFunc[analysis]

            return funcs if type(funcs) == list else [funcs]
        except KeyError:
            return None

    def analyze(self, repertoire, analyses=Ac.ALL):
        if analyses == Ac.ALL:
            analyses = self.analysisToFunc.keys()

        # For each analysis type specified in the list,
        for analysis in analyses:
            if self.VERBOSE and not self.parallel:
                sys.stdout.write(Ac.analysisToDesc[analysis] + " ... ")

            # Get a list of function names
            functions = self.getFuncsFor(analysis)

            # For each function name,
            for function in functions:
                # Call the function for the given repertoire
                function(repertoire)

        if self.VERBOSE and not self.parallel:
            sys.stdout.write("Done.")

    def getLandscapeObj(self, giant, repertoire):
        # Get the landscape object
        lscape = Landscape(giant, self.netBuilder, self.inDataDict[repertoire],
                           self.deltaDict[repertoire], self.bitManip)

        return lscape

    def peaks(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Get the landscape object
        lscape = self.getLandscapeObj(giant, repertoire)

        # Get peaks
        peaks = lscape.getPeaks(recompute=True)

        # Set the computed values as a network attribute
        giant["Number_of_peaks"] = len(peaks)

        # Store a dict - {key=peakId, value=[sequences in the peak]}
        giant["Peaks"] = {
            peakId: peaks[peakId]["sequences"]
            for peakId in peaks.keys()
        }

        # For each vertex, add distance to summit as a vertex attribute
        lscape.populateDistsToSummit()

    def paths(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Get the landscape object
        lscape = self.getLandscapeObj(giant, repertoire)

        # Get paths
        lscape.getAccessiblePaths(0)

        # Set the computed value as a network attribute
        giant["Summit"] = lscape.pathAnalyzer.getSummitId()

        # Vertex level attributes
        allPathsToPeak = lscape.pathAnalyzer.getAllPathsToPeak()
        giant.vs["pathsToSummit"] = [
            allPathsToPeak[i]
            for i in range(len(allPathsToPeak))
        ]

        # Add the count for each vertex as a vertex level attribute in giant
        giant.vs["Accessible_paths_through"] = lscape.pathAnalyzer.getPathsThruVtxs()

    def paths_ratios(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # 'PathAnalyzer' object
        path_analyzer = PathAnalyzer(giant, self.netBuilder,
                                     self.deltaDict[repertoire])

        # Run shortest paths calculation for all paths; regardless of path length.
        # This sets the 'max_path_length' value.
        path_analyzer.getAccessiblePaths()

        # Length of the longest path in the network
        max_path_length = path_analyzer.max_path_length

        # Compute the ratio of accessible paths for all path lengths in range:
        # [2, max_path_length].
        # Set dict {path_length : ratio}
        giant["Ratio_of_accessible_mutational_paths"] = {
            i: path_analyzer.getAccessiblePaths(i)
            for i in xrange(2, max_path_length + 1)
        }

    def epistasis(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Get the landscape object
        lscape = self.getLandscapeObj(giant, repertoire)

        # Get epistasis
        epistasis = lscape.getEpistasis()

        # Get the vertex to squares dict
        vtxToSqrs, squares = lscape.epiAnalyzer.getVertexToSquaresDict()

        # Set the computed values as network attributes
        giant["Squares_list"] = json.dumps(squares)
        giant["SqrEpi_list"] = json.dumps(lscape.epiAnalyzer.getSqrEpi())
        giant["Number_of_squares"] = len(lscape.epiAnalyzer.squares)
        giant["Magnitude_epistasis"] = epistasis[Epi.MAGNITUDE]
        giant["Simple_sign_epistasis"] = epistasis[Epi.SIGN]
        giant["Reciprocal_sign_epistasis"] = epistasis[Epi.RECIPROCAL_SIGN]

        # Vertex level attributes

        # Can't traverse an empty dict
        if len(vtxToSqrs) > 0:
            # Since vtxToSqrs is a dict, it is important to convert it into an ordered
            # list for correct mapping to vertices
            giant.vs["VtxToSqrs"] = [vtxToSqrs[i] for i in range(len(vtxToSqrs))]
        else:
            # Empty list. This makes it possible for the user to distinguish between
            # no squares associated with a vertex, and epistasis function was not
            # called at all
            # Note: Has to be a nested empty list, since otherwise igraph throws an
            # exception
            giant.vs["VtxToSqrs"] = [[]]

    def robustness(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Construct a RobustnessAnalyzer object
        robAnalyzer = RobustnessAnalyzer(giant, self.netBuilder)

        # Compute repertoire robustness and set it as a network
        # attribute
        giant["Robustness"] = robAnalyzer.getAvgRobustness()

        # Set robustness values for all vertices, i.e., sequences
        giant.vs["Robustness"] = robAnalyzer.getRobustnessAll()

    # Data structure initializations that need only be done once for
    # evolvability analysis of all repertoires.
    def init_evolvability(self):
        self.seqToRepDict_evo = EvolvabilityAnalyzer.updateSeqToRepDict(
            self.seqToRepDict, self.repToGiantDict)

        if self.isDoubleStranded:
            self.rcToSeqDict = EvolvabilityAnalyzer.buildRcToSeqDict(
                self.seqToRepDict_evo, self.bitManip)

        self.bitsToSeqDict = EvolvabilityAnalyzer.buildBitsToSeqDict(
            self.seqToRepDict_evo, self.rcToSeqDict, self.bitManip,
            self.isDoubleStranded)

    def evolvability(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Construct a EvolvabilityAnalyzer object
        evoAnalyzer = EvolvabilityAnalyzer(giant,
                                           self.inDataDict,
                                           self.seqToRepDict_evo,
                                           self.repToGiantDict,
                                           self.rcToSeqDict,
                                           self.bitsToSeqDict,
                                           self.netBuilder,
                                           self.isDoubleStranded)

        # Compute repertoire evolvability and set it as a network
        # attribute
        repertoireEvo, targetRepertoires = evoAnalyzer.getReportoireEvo()
        giant["Evolvability"] = repertoireEvo

        # Stringify the list, since pythons lists cannot be written to GML.
        giant["Evolvability_targets"] = json.dumps(targetRepertoires)

        # Set evolvability values for all vertices, i.e., sequences
        evoTuples = evoAnalyzer.getEvoAll()

        evoScores = [evoTuples[i][0] for i in range(len(evoTuples))]
        evoTargets = [evoTuples[i][1] for i in range(len(evoTuples))]

        giant.vs["Evolvability"] = evoScores
        giant.vs["Evolves_to_genotypes_in"] = [evoTargets[i].keys()
                                               for i in range(len(evoTargets))]
        giant.vs["Evolvability_targets"] = evoTargets

    def accessibility(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Create an AccessibilityAnalyzer object
        accAnalyzer = AccessibilityAnalyzer(repertoire, giant,
                                            self.repToGiantDict,
                                            self.inDataDict,
                                            self.netBuilder,
                                            self.bitManip,
                                            self.isDoubleStranded)

        # Compute repertoire accessibility and set it as a network attribute
        giant["Accessibility"] = accAnalyzer.getAccessibility()

    def neighborAbundance(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Create an AccessibilityAnalyzer object
        accAnalyzer = AccessibilityAnalyzer(repertoire, giant,
                                            self.repToGiantDict,
                                            self.inDataDict,
                                            self.netBuilder,
                                            self.bitManip,
                                            self.isDoubleStranded)

        # Compute repertoire neighborhood abundance and set it as a network
        # attribute
        giant["Neighbor_abundance"] = accAnalyzer.getNeighborAbundance()

    def phenotypicDiversity(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Create an AccessibilityAnalyzer object
        accAnalyzer = AccessibilityAnalyzer(repertoire, giant,
                                            self.repToGiantDict,
                                            self.inDataDict,
                                            self.netBuilder,
                                            self.bitManip,
                                            self.isDoubleStranded)

        # Compute phenotypic diversity and set it as a network attribute
        giant["Diversity_index"] = accAnalyzer.getPhenotypicDivesity()

    def structure(self, repertoire):
        # Get the genotype network for the repertoire
        network = self.caller.genotype_network(repertoire)

        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Create the structure analyzer object
        structAnalyzer = StructureAnalyzer(network, self.netBuilder)

        # Compute and set network/giant level properties
        network["Genotype_network_sizes"] = str(structAnalyzer.getComponentSizes())
        network["Number_of_genotype_networks"] = structAnalyzer.getNumComponents()
        network["Size_of_dominant_genotype_network"] = structAnalyzer.getDominantSize()
        network["Proportional_size_of_dominant_genotype_network"] = structAnalyzer.getPercentDominantSize()
        giant["Edge_density"] = structAnalyzer.getEdgeDensity()
        giant["Diameter"] = structAnalyzer.getDiameter()
        giant["Average_clustering_coefficient_of_dominant_genotype_network"] = structAnalyzer.getAvgClstrCoeff()
        giant["Assortativity"] = structAnalyzer.getAssortativity()

        # The list of vertex Ids needs to be stringified, since otherwise these
        # cannot be written to GML.
        giant["diameterPath_list"] = json.dumps(structAnalyzer.getDiameterPath())

        # Compute and set vertex level properties
        giant.vs["Coreness"] = structAnalyzer.getCoreness()
        giant.vs["Clustering_coefficient"] = structAnalyzer.getClusteringCoefficients()

    # The parameter 'r' is just a place holder, and is needed in the
    # signature just so that it can be called anonymously from
    # analyze(). See analyze() for further details.
    def overlap(self, r=None):
        # Overlap analysis is not allowed during parallel processing,
        # since it needs to be performed only once
        if self.parallel:
            return

        # If overlap has been computed for any of the repertoires already,
        # it does not need to be computed again.
        if not self.overlapMatrix:
            # Create the overlap analyzer
            overlapAnalyzer = OverlapAnalyzer(self.repToGiantDict,
                                              self.caller.genotype_sets())

            # Compute overlap data
            self.overlapMatrix, repertoires, overlapDict = overlapAnalyzer.getOverlapData()

            # If the overlap dict was populated
            if overlapDict:
                # Use the overlap dict to populate vertex level attributes in all giants

                # For each repertoire,
                for repertoire in self.caller.genotype_sets():
                    # Get giant
                    giant = self.repToGiantDict[repertoire]

                    # List of all unique repertoires that overlap with this giant
                    giant["Overlapping_genotype_sets"] = set()

                    # Get the sequence dict for this repertoire
                    seqDict = overlapDict[repertoire]

                    # For each sequence in seq dict,
                    for sequence in seqDict.keys():
                        # Get the corresponding vertex id from the network
                        vertex = self.netBuilder.getVertex(sequence, giant)

                        # List of repertoires that contain the sequence
                        overlapping_seqs = seqDict[sequence]

                        # Add the list of targets as an attribute for this
                        # vertex in giant
                        giant.vs[vertex.index]["Overlaps_with_genotypes_in"] = \
                            overlapping_seqs

                        # Add the overlapping repertoires to the set of all
                        # repertoires that overlap with any sequence in this
                        # repertoire.
                        giant["Overlapping_genotype_sets"] |= set(overlapping_seqs)

                    # Convert the set to list for easy output file writing
                    giant["Overlapping_genotype_sets"] = \
                        list(giant["Overlapping_genotype_sets"])

                    # Calculate the ratio of No. of overlapping repertoires to
                    # the total No. of other repertoires
                    try:
                        ratio = float(len(giant["Overlapping_genotype_sets"])) / \
                                (float(len(self.caller.genotype_sets())) - 1)
                    except ZeroDivisionError:
                        ratio = 0

                    giant["Ratio_of_overlapping_genotype_sets"] = ratio

            # If overlap matrix was populated,
            if self.overlapMatrix:
                # Write matrix to file
                Writer.writeOverlapToFile(self.overlapMatrix, repertoires,
                                          self.caller.cmdArgs.outPath)
