"""
    analysis_handler
    ~~~~~~~~~~~~~~~~

    Exposes wrapper functions, one per analysis type. Serves as a collective
    interface to handler classes for all analysis types.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import json  # For proper list stringification

from genonets.writer import Writer
from genonets.filters import WriterFilter
from genonets.peak_functions import PeakAnalyzer
from genonets.path_functions import PathAnalyzer
from genonets.overlap_functions import OverlapAnalyzer
from genonets.structure_functions import StructureAnalyzer
from genonets.epistasis_functions import EpistasisAnalyzer
from genonets.epistasis_functions import AlternateNetworkData
from genonets.robustness_functions import RobustnessAnalyzer
from genonets.evolvability_functions import EvolvabilityAnalyzer
from genonets.accessibility_functions import AccessibilityAnalyzer
from genonets.constants import AnalysisConstants as Ac
from genonets.constants import EpistasisConstants as Epi


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

        self._alternate_bit_manip = self.caller.alternate_bit_manip
        self._alternate_net_builder = self.caller.alternate_net_builder
        self._rep_to_alternate_net_dict = self.caller.rep_to_alternate_net_dict
        self._rep_to_alternate_giant_dict = self.caller.rep_to_alternate_giant_dict

        # Dictionary to store 'analysis type' to 'function'
        # mapping
        self.analysisToFunc = {
            Ac.PEAKS: self.peaks,
            Ac.PATHS: self.paths,
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
        # reverse complements should be used in various analyses
        self.isDoubleStranded = self.caller.cmd_args.use_reverse_complements

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
                print(
                    f'{Ac.analysisToDesc[analysis]}', end=' ... \n', flush=True
                )

            # Get a list of function names
            functions = self.getFuncsFor(analysis)

            # For each function name,
            for function in functions:
                # Call the function for the given repertoire
                function(repertoire)

        if self.VERBOSE and not self.parallel:
            print('Done.', flush=True)

    def peaks(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Create the peak analyzer
        analyzer = PeakAnalyzer(
            giant, self.netBuilder, self.deltaDict[repertoire], self.VERBOSE)

        # Get peaks
        peaks = analyzer.get_all_peaks(recompute=True)

        # Convert sets to lists for a simpler output format
        for peak_id in peaks:
            peaks[peak_id] = list(peaks[peak_id])

        # Set the computed values as network attributes
        giant['Peaks'] = peaks
        giant['Number_of_peaks'] = len(peaks)

    def paths(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        if giant.vcount() < 2:
            print(f'Warning: Cancelling paths analysis. Paths can only be '
                  f'computed if there are at least two vertices in the '
                  f'network.')

            return

        # 'PathAnalyzer' object
        path_analyzer = PathAnalyzer(
            giant,
            self.netBuilder,
            self.deltaDict[repertoire],
            self.VERBOSE,
            self.inDataDict[repertoire]
        )

        result_shortest_paths = path_analyzer.compute_shortest_paths()
        result_indirect_paths = path_analyzer.compute_indirect_paths()

        giant["Summit"] = path_analyzer.summit_sequences

        giant["Ratio_of_accessible_mutational_paths"] = result_shortest_paths[
            "Ratio_of_accessible_mutational_paths"]
        giant.vs["Accessible_paths_from"] = result_shortest_paths[
            "Accessible_paths_from"]
        giant.vs["Shortest_path_length"] = result_shortest_paths[
            "Shortest_path_length"]

        giant.vs["Shortest_accessible_path_length"] = result_indirect_paths

    def epistasis(self, repertoire):
        # Dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        alternate_net_data = None

        if self.caller.letter_to_neighbors:
            self.caller.create_alternate_network(repertoire)

            self._rep_to_alternate_net_dict = \
                self.caller.rep_to_alternate_net_dict
            self._rep_to_alternate_giant_dict = \
                self.caller.rep_to_alternate_giant_dict

            alternate_net_data = AlternateNetworkData(
                self._rep_to_alternate_giant_dict[repertoire],
                self._alternate_bit_manip,
                self._alternate_net_builder
            )

        # Epistasis analyzer
        epi_analyzer = EpistasisAnalyzer(
            network=giant,
            net_utils=self.netBuilder,
            genotype_to_score_map=self.inDataDict[repertoire],
            delta=self.deltaDict[repertoire],
            alternate_net_data=alternate_net_data,
            verbose=self.VERBOSE,
            save_squares=self.caller.cmd_args.save_squares,
            out_dir=self.caller.cmd_args.output_path
        )

        # Compute epistasis
        epistasis = epi_analyzer.compute_epistasis_ratios()

        # Set the computed values as network attributes
        giant['Number_of_squares'] = epi_analyzer.num_squares
        giant['Magnitude_epistasis'] = epistasis[Epi.MAGNITUDE]
        giant['Simple_sign_epistasis'] = epistasis[Epi.SIGN]
        giant['Reciprocal_sign_epistasis'] = epistasis[Epi.RECIPROCAL_SIGN]

    def robustness(self, repertoire):
        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Construct a RobustnessAnalyzer object
        robAnalyzer = RobustnessAnalyzer(
            giant,
            self.netBuilder,
            self.isDoubleStranded,
            self.VERBOSE
        )

        # Compute repertoire robustness and set it as a network
        # attribute
        giant["Robustness"] = robAnalyzer.getAvgRobustness()

        # Set robustness values for all vertices, i.e., sequences
        giant.vs["Robustness"] = robAnalyzer.getRobustnessAll()

    # Data structure initializations that need only be done once for
    # evolvability analysis of all repertoires.
    def init_evolvability(self):
        if self.caller.cmd_args.use_all_components:
            # Use the entire entire network
            self.seqToRepDict_evo = self.seqToRepDict
        else:
            # Use only the dominant network, i.e, the giant component
            self.seqToRepDict_evo = EvolvabilityAnalyzer.updateSeqToRepDict(
                self.seqToRepDict, self.repToGiantDict
            )

        if self.isDoubleStranded:
            self.rcToSeqDict = EvolvabilityAnalyzer.buildRcToSeqDict(
                self.seqToRepDict_evo, self.bitManip
            )

        self.bitsToSeqDict = EvolvabilityAnalyzer.buildBitsToSeqDict(
            self.seqToRepDict_evo,
            self.rcToSeqDict,
            self.bitManip,
            self.isDoubleStranded
        )

    def evolvability(self, repertoire):
        # Evolvability analysis only makes sense if there are at least two
        # repertoires
        if len(self.repToNetDict) < 2:
            print('Cancelling evolvability analysis, as it requires at '
                  'least two genotype sets')
            return

        if self.caller.cmd_args.use_all_components:
            # Use the entire entire network
            network = self.caller.genotype_network(repertoire)
        else:
            # Use only the dominant network, i.e, the giant component
            network = self.caller.dominant_network(repertoire)

        # Construct a EvolvabilityAnalyzer object
        evolvability_analyzer = EvolvabilityAnalyzer(
            network,
            self.inDataDict,
            self.seqToRepDict_evo,
            self.repToGiantDict,
            self.rcToSeqDict,
            self.bitsToSeqDict,
            self.netBuilder,
            self.isDoubleStranded,
            self.VERBOSE
        )

        # Compute repertoire evolvability and set it as a network
        # attribute
        repertoire_evolvability, target_repertoires = \
            evolvability_analyzer.getReportoireEvo()
        network["Evolvability"] = repertoire_evolvability

        # Stringify the list, since pythons lists cannot be written to GML.
        network["Evolvability_targets"] = json.dumps(sorted(target_repertoires))

        # Set evolvability values for all vertices, i.e., sequences
        evolvability_tuples = evolvability_analyzer.getEvoAll()

        evolvability_scores = [
            evolvability_tuples[i][0] for i in range(len(evolvability_tuples))
        ]
        evolvability_targets = [
            evolvability_tuples[i][1] for i in range(len(evolvability_tuples))
        ]

        network.vs["Evolvability"] = evolvability_scores
        network.vs["Evolves_to_genotypes_in"] = [
            sorted([*evolvability_targets[i]])
            for i in range(len(evolvability_targets))
        ]
        network.vs["Evolvability_targets"] = evolvability_targets

        network['Interface_edges'] = evolvability_analyzer.compute_external_mutation_types(
            target_repertoires=target_repertoires,
            targets_per_genotype=evolvability_targets
        )

    def accessibility(self, repertoire):
        # Accessibility analysis only makes sense if there are at least
        # two repertoires
        if len(self.repToGiantDict) < 2:
            print('Cancelling accessibility analysis, as it requires at '
                  'least two genotype sets')
            return

        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Create an AccessibilityAnalyzer object
        accAnalyzer = AccessibilityAnalyzer(
            repertoire,
            giant,
            self.repToGiantDict,
            self.inDataDict,
            self.netBuilder,
            self.bitManip,
            self.isDoubleStranded,
            self.VERBOSE
        )

        # Compute repertoire accessibility and set it as a network attribute
        giant["Accessibility"] = accAnalyzer.getAccessibility()

    def neighborAbundance(self, repertoire):
        # Neighbor abundance analysis only makes sense if there are at least
        # two repertoires
        if len(self.repToGiantDict) < 2:
            print('Cancelling neighbor abundance analysis, as it requires at '
                  'least two genotype sets')
            return

        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Create an AccessibilityAnalyzer object
        accAnalyzer = AccessibilityAnalyzer(repertoire, giant,
                                            self.repToGiantDict,
                                            self.inDataDict,
                                            self.netBuilder,
                                            self.bitManip,
                                            self.isDoubleStranded,
                                            self.VERBOSE)

        # Compute repertoire neighborhood abundance and set it as a network
        # attribute
        giant["Neighbor_abundance"] = accAnalyzer.getNeighborAbundance()

    def phenotypicDiversity(self, repertoire):
        # Phenotypic diversity analysis only makes sense if there are at least
        # two repertoires
        if len(self.repToGiantDict) < 2:
            print('Cancelling phenotypic diversity analysis, as it requires at '
                  'least two genotype sets')
            return

        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Create an AccessibilityAnalyzer object
        accAnalyzer = AccessibilityAnalyzer(repertoire, giant,
                                            self.repToGiantDict,
                                            self.inDataDict,
                                            self.netBuilder,
                                            self.bitManip,
                                            self.isDoubleStranded,
                                            self.VERBOSE)

        # Compute phenotypic diversity and set it as a network attribute
        giant["Diversity_index"] = accAnalyzer.getPhenotypicDivesity()

    def structure(self, repertoire):
        # Get the genotype network for the repertoire
        network = self.caller.genotype_network(repertoire)

        # Get the dominant genotype network for the repertoire
        giant = self.caller.dominant_network(repertoire)

        # Create the structure analyzer object
        structAnalyzer = StructureAnalyzer(network, self.netBuilder, self.VERBOSE)

        # Compute and set network/giant level properties
        network["Genotype_network_sizes"] = str(structAnalyzer.getComponentSizes())
        network["Number_of_genotype_networks"] = structAnalyzer.getNumComponents()
        network["Size_of_dominant_genotype_network"] = structAnalyzer.getDominantSize()
        network["Proportional_size_of_dominant_genotype_network"] = structAnalyzer.getPercentDominantSize()
        giant["Edge_density"] = structAnalyzer.getEdgeDensity()
        # if self.VERBOSE:
        #     print('Retrieving diameter ...')
        # giant["Diameter"] = structAnalyzer.getDiameter()
        giant["Average_clustering_coefficient_of_dominant_genotype_network"] = structAnalyzer.getAvgClstrCoeff()
        giant["Assortativity"] = structAnalyzer.getAssortativity()

        # # The list of vertex Ids needs to be stringified, since otherwise these
        # # cannot be written to GML.
        # giant["diameterPath_list"] = json.dumps(structAnalyzer.getDiameterPath())

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
                                              self.caller.genotype_sets(),
                                              self.bitManip,
                                              self.isDoubleStranded,
                                              WriterFilter.genotype_set_to_order)

            # Compute overlap data. Note: The list of genotype sets is returned from
            # the function to make sure the order used inside the function is the one
            # used for further calculations here.
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
                        sorted(list(giant["Overlapping_genotype_sets"]))

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
                Writer.writeOverlapToFile(
                    self.overlapMatrix,
                    repertoires,
                    self.caller.cmd_args.output_path
                )
