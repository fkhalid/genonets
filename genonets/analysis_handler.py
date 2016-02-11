"""
    analysis_handler
    ~~~~~~~~~~~~~~~~

    Exposes wrapper functions, one per analysis type. Serves as a collective interface to handler classes for all analysis types.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import json  # For proper list stringification

from genonets_writer import Writer
from landscape_functions import Landscape
from overlap_functions import OverlapAnalyzer
from structure_functions import StructureAnalyzer
from robustness_functions import RobustnessAnalyzer
from evolvability_functions import EvolvabilityAnalyzer
from accessibility_functions import AccessibilityAnalyzer
from genonets_constants import AnalysisConstants as ac
from genonets_constants import EpistasisConstants as epi


class AnalysisHandler:
    # Constructor
    def __init__(self, caller, parallel=False):
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

        # Dict {bitseq : seq}
        self.bitsToSeqDict = self.caller.bitsToSeqDict

        # Dictionary to store 'analysis type' to 'function'
        # mapping
        self.analysisToFunc = {  # ac.LANDSCAPE : self.landscape, \
            ac.PEAKS: self.peaks,
            ac.PATHS: self.paths,
            ac.EPISTASIS: self.epistasis,
            ac.ROBUSTNESS: self.robustness,
            ac.EVOLVABILITY: self.evolvability,
            ac.ACCESSIBILITY: self.accessibility,
            ac.NEIGHBOR_ABUNDANCE: self.neighborAbundance,
            ac.PHENOTYPIC_DIVERSITY: self.phenotypicDiversity,
            ac.STRUCTURE: self.structure,
            ac.OVERLAP: self.overlap
        }

        # Refernce to the overlap matrix
        self.overlapMatrix = None

    def getFuncsFor(self, analysis):
        try:
            funcs = self.analysisToFunc[analysis]

            return funcs if type(funcs) == list else [funcs]
        except KeyError:
            return None

    def analyze(self, repertoire, analyses=ac.ALL):
        if analyses == ac.ALL:
            analyses = self.analysisToFunc.keys()

        # For each analysis type specified in the list,
        for analysis in analyses:
            # Get a list of function names
            functions = self.getFuncsFor(analysis)

            # For each function name,
            for function in functions:
                # Call the function for the given repertoire
                function(repertoire)

    def getLandscapeObj(self, giant, repertoire):
        # Get the landscape object
        lscape = Landscape(giant, self.netBuilder, self.inDataDict[repertoire],
                           self.deltaDict[repertoire], self.bitManip)

        return lscape

    def peaks(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.getDominantNetFor(repertoire)

        # Get the landscape object
        lscape = self.getLandscapeObj(giant, repertoire)

        # Get peaks
        peaks = lscape.getPeaks(recompute=True)

        # Set the computed values as a network attribute
        giant["Number_of_peaks"] = len(peaks)

        # Store a dict - {key=peakId, value=[sequences in the peak]}
        giant["Peaks"] = {peakId: peaks[peakId]["sequences"] \
                          for peakId in peaks.keys()}

        # For each vertex, add distance to summit as a vertex attribute
        lscape.populateDistsToSummit()

    def paths(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.getDominantNetFor(repertoire)

        # Get the landscape object
        lscape = self.getLandscapeObj(giant, repertoire)

        # Get paths
        lscape.getAccessiblePaths(0)

        # Set the computed value as a network attribute
        giant["Summit"] = lscape.pathAnalyzer.getSummitId()

        # Vertex level attributes
        allPathsToPeak = lscape.pathAnalyzer.getAllPathsToPeak()
        giant.vs["pathsToSummit"] = [allPathsToPeak[i] \
                                     for i in range(len(allPathsToPeak))]

        # Add the count for each vertex as a vertex level attribute in giant
        giant.vs["Accessible_paths_through"] = lscape.pathAnalyzer.getPathsThruVtxs()

    def epistasis(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.getDominantNetFor(repertoire)

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
        giant["Magnitude_epistasis"] = epistasis[epi.MAGNITUDE]
        giant["Simple_sign_epistasis"] = epistasis[epi.SIGN]
        giant["Reciprocal_sign_epistasis"] = epistasis[epi.RECIPROCAL_SIGN]

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
        giant = self.caller.getDominantNetFor(repertoire)

        # Construct a RobustnessAnalyzer object
        robAnalyzer = RobustnessAnalyzer(giant, self.netBuilder)

        # Compute repertoire robustness and set it as a network
        # attribute
        giant["Robustness"] = robAnalyzer.getAvgRobustness()

        # Set robustness values for all vertices, i.e., sequences
        giant.vs["Robustness"] = robAnalyzer.getRobustnessAll()

    def evolvability(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.getDominantNetFor(repertoire)

        # Construct a EvolvabilityAnalyzer object
        evoAnalyzer = EvolvabilityAnalyzer(giant, self.inDataDict,
                                           self.seqToRepDict, self.repToGiantDict, self.netBuilder,
                                           self.bitsToSeqDict)

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
        giant = self.caller.getDominantNetFor(repertoire)

        # Create an AccessibilityAnalyzer object
        accAnalyzer = AccessibilityAnalyzer(repertoire, giant, self.repToGiantDict,
                                            self.inDataDict, self.netBuilder, self.bitManip)

        # Compute repertoire accessibility and set it as a network attribute
        giant["Accessibility"] = accAnalyzer.getAccessibility()

    def neighborAbundance(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.getDominantNetFor(repertoire)

        # Create an AccessibilityAnalyzer object
        accAnalyzer = AccessibilityAnalyzer(repertoire, giant, self.repToGiantDict,
                                            self.inDataDict, self.netBuilder, self.bitManip)

        # Compute repertoire neighborhood abundance and set it as a network
        # attribute
        giant["Neighbor_abundance"] = accAnalyzer.getNeighborAbundance()

    def phenotypicDiversity(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.getDominantNetFor(repertoire)

        # Create an AccessibilityAnalyzer object
        accAnalyzer = AccessibilityAnalyzer(repertoire, giant, self.repToGiantDict,
                                            self.inDataDict, self.netBuilder, self.bitManip)

        # Compute phenotypic diversity and set it as a network attribute
        giant["Phenotypic_diversity"] = accAnalyzer.getPhenotypicDivesity()

    def structure(self, repertoire):
        # Get the genotype network for the repertoire
        network = self.caller.getNetworkFor(repertoire)

        # Get the dominant genotype network for the repertoire
        giant = self.caller.getDominantNetFor(repertoire)

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

    def overlap(self, r=None):
        # Overlap analysis is not allowed during parallel processing,
        # since it needs to be performed only once
        if self.parallel == True:
            return

        # If overlap has been computed for any of the repertoires already,
        # it does not need to be computed again.
        if not self.overlapMatrix:
            # Creat the overlap analyzer
            overlapAnalyzer = OverlapAnalyzer(self.repToGiantDict,
                                              self.caller.getRepertoires())

            # Compute overlap data
            self.overlapMatrix, repertoires, overlapDict = overlapAnalyzer.getOverlapData()

            # If the overlap dict was populated
            if overlapDict:
                # Use the overlap dict to populate vertex level attributes in all giants

                # For each repertoire,
                for repertoire in self.caller.getRepertoires():
                    # Get giant
                    giant = self.repToGiantDict[repertoire]

                    # Get the sequence dict for this repertoire
                    seqDict = overlapDict[repertoire]

                    # For each sequence in seq dict,
                    for sequence in seqDict.keys():
                        # Get the corresponding vertex id from the network
                        vertex = self.netBuilder.getVertex(sequence, giant)

                        # Add the list of targets as an attribute for this
                        # vertex in giant
                        giant.vs[vertex.index]["Overlaps_with_genotypes_in"] = \
                            seqDict[sequence]

            # If overalp matrix was populated,
            if self.overlapMatrix:
                # Write matrix to file
                Writer.writeOverlapToFile(self.overlapMatrix, repertoires,
                                          self.caller.cmdArgs.outPath)
