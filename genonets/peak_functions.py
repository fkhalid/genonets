"""
    This module implements the functions required for identification of peaks
    in a network (when viewed as a fitness landscape).

    Note: The term 'peak' is used in this module to refer to both peaks and
          plateaus.

    The algorithm:
        1. Sort all genotypes in descending order of score.
        2. For each genotype,
        3.      Create a list of all genotypes that are in the same score band
                as the focal genotype. Exclude genotypes that have already been
                processed.
        4.      Create the neutral zone, i.e., a list of all genotypes that are
                in the score band and connected. The focal genotype is also
                included in the neutral zone.
        5.      If any genotype in the neutral zone has a 1-neighbor outside
                the neutral zone that has already been processed (either as a
                peak member of as a valley member), the entire neutral zone is
                a non-peak.
        6.      Else, create a new peak such that all genotypes in the neutral
                zone are members of the peak.
        7.      Mark all genotypes in the neutral zone as processed.

"""

from genonets_utils import Utils


class PeakAnalyzer:
    """
    Encapsulates the computations required for identifying peaks in the given
    network. It is assumed that the given network is a single connected
    component; typically the dominant component of a genotype network.

    """

    def __init__(self, network, net_utils, delta):
        """
        Initializes the object.

        :param network: (igraph.Graph) The network within which to identify
                peaks. The network must be connected, i.e., it must not consist
                of more than one connected component.
        :param net_utils: (genonets.graph_utils.NetworkBuilder) Reference to
                the object providing the required network manipulation
                functions.
        :param delta: (float) The experimental noise value to consider when
                evaluating each genotype for optimality.

        """

        self._delta = delta
        self._network = network
        self._net_utils = net_utils

        # Reference to the BitSeqManipulator in use
        self._bit_manip = net_utils.bitManip

        # Map { genotype: {all 1-neighbors} }
        self._neighbor_map = {
            s: set(self._net_utils.getNeighborSequences(s, self._network))
            for s in network.vs['sequences']
        }

        # Set of elements that have already been processed, either individually
        # or as part of a neutral zone for another genotype
        self._processed_genotypes = set()

        # Map { peak_id: {genotypes} }. To be populated with all computed peaks.
        self._peaks = None

    def get_peak_with_genotype(self, genotype, recompute=False):
        """
        Returns the peak object of which the given genotype is a member.

        :param genotype: (str) The genotype to search in the peaks.
        :param recompute: (bool) 'True' if peaks should be computed again even
                if the corresponding analyzer has been used to compute peaks
                already. 'False' otherwise.

        :return: (dict) of the form {peak_id: {all genotypes in the peak}} for
                the peak in which the given genotype was found.
                (None) if the genotype could not be found in any peak. This
                would be the case if the genotype lies in a valley.

        """

        # Get all the peaks
        peaks = self.get_all_peaks(recompute)

        for peak_id in peaks:
            # Get the set of genotypes that belong to this peak
            peak_members = peaks[peak_id]

            # If the given genotype is a member,
            if genotype in peak_members:
                return {peak_id: peak_members}

        return None

    def get_all_peaks(self, recompute=False):
        """
        Returns all peaks in the network with which the corresponding analyzer
        object was initialized.

        If the corresponding analyzer object has been used to compute peaks
        already, and 'recompute' is 'False', the already computed peaks object
        is returned. Otherwise, peaks are computed first and then returned.

        :param recompute: (bool) 'True' if peaks should be computed again even
                if the corresponding analyzer has been used to compute peaks
                already. 'False' otherwise.

        :return: (dict) of the form {peak_id: {all genotypes in the peak}}.

        """

        if not self._peaks or recompute:
            # Structured array of tuples (genotype, score), sorted in
            # descending order of scores
            genotype_score_pairs = Utils.getSortedSeqEscArr(
                self._network,
                self._bit_manip.seqLength,
                sortOrder='descending'
            )

            # Identify peaks
            self._peaks = self._build_peaks(genotype_score_pairs)

        return self._peaks

    def _build_peaks(self, elements):
        """
        Assuming the given array of (genotype, score) tuples is already sorted
        in descending order of scores, this function orchestrates the algorithm
        for determining whether a given genotype lies in a peak or a valley.
        That is, this function implements steps 2 - 7 of the overall algorithm
        specified in the documentation of this module.

        :param elements: (list of tuples), where each tuple is of the form
                    (genotype, score). The list is assumed to be in descending
                    order of scores.

        :return: (dict) { key=peak_id: value={all genotypes in the peak} }

        """

        peaks = {}

        for i in xrange(len(elements)):
            # If the current element has already been processed,
            if elements[i]['sequence'] in self._processed_genotypes:
                # Move on to the next element
                continue

            # List of elements for which the score lies within the score band
            # defined by the score of the current focal element
            score_band = self._build_score_band(i, elements)

            # Neutral zone, i.e., all elements within the score band that are
            # connected
            neutral_zone = self._build_neutral_zone(elements[i], score_band)

            # Append the focal element itself to the neutral zone.
            # Note: Insert at index 0 would be nicer, but the cost would O(n).
            neutral_zone.append(elements[i])

            # If any genotype in the neutral zone has a 1-neighbor outside
            # the neutral zone that has already been processed (either as a
            # peak member of as a valley member),
            if self._neighbors_processed_element(neutral_zone):
                # The entire neutral zone is a non-peak. There's nothing more
                # to be done.
                pass
            else:
                # Create a new peak, and add all members of the neutral zone
                # to this peak
                peak_id = len(peaks)
                peaks[peak_id] = self._create_peak(neutral_zone)

            # Mark all elements in the neutral zone as processed, i.e., there
            # is no need to process each one individually
            self._processed_genotypes.update({
                e['sequence'] for e in neutral_zone
            })

        return peaks

    def _build_score_band(self, index_of_focal, elements):
        """
        Returns a list of elements that fall within the noise determined score 
        band for the given focal element (excluding elements that have already
        been processed), where the focal element is identified by the given
        index.

        :param index_of_focal: (int) Index of the focal element in the given
                list of all elements.
        :param elements: (numpy.array) Array of (genotype, score) tuples
                corresponding to all genotypes in the network in consideration.
        
        :return: (list) Elements for which the score lies within the band of
                the focal element, and which have not been processed already.
        
        """

        score_band = []
        focal_element = elements[index_of_focal]

        # Starting from the element succeeding the focal element,
        for i in xrange(index_of_focal + 1, len(elements)):
            # If the score lies within the band,
            if elements[i]['escore'] >= focal_element['escore'] - self._delta:
                # If the current element has already been processed, it means
                # that it is either part of an existing peak, or an existing
                # non-peak. It can only be indirectly connected to the current
                # score band via an existing peak, or via an existing non-peak.
                # Either way, we can safely ignore it, as other checks should
                # be sufficient to accurately categorize the current score band.
                if elements[i]['sequence'] in self._processed_genotypes:
                    pass
                else:
                    score_band.append(elements[i])
            else:
                # Once we've reached an element which lies in a lower band, we
                # know that the remaining elements will also lie in lower
                # bands. Therefore, we don't need to look any further.
                break

        return score_band

    def _build_neutral_zone(self, focal_element, score_band):
        """
        Returns the neutral zone for the given focal element.

        The neutral zone consists of all elements that satisfy the following
        conditions:
            1) The element lies within the score band for the given focal 
               element.
            2) The element is either a 1-neighbor of the focal element, or it
               is indirectly connected to the focal element, such that all
               elements that constitute the path between this element and the
               focal element lie within the score band of the focal element.
        
        :param focal_element: (tuple) of form (genotype, score) corresponding
                to the element for which the score band was constructed.
        :param score_band: (list) of (genotype, score) tuples, where each
                tuple represents a genotype for which the score lies in the
                score band for the focal element.
        
        :return: (list) The neutral zone.
        
        """

        # 1-neighbors of the focal element
        neighbors_of_focal = self._get_neighbors_of(focal_element['sequence'])

        # Elements in the given score band that are also in the neighborhood of
        # the focal element. These are by definition already part of the neutral
        # zone.
        neutral_zone = [
            item for item in score_band
            if item['sequence'] in neighbors_of_focal
        ]

        # Elements in the given score band that are not 1-neighbors of the
        # focal element
        band_non_neighbors = [e for e in score_band if e not in neutral_zone]

        # Perform a breadth-first search of the neutral zone to find those
        # elements in band-non-neighbors that might be indirectly connected
        # (via other members of the neutral zone) to the focal element.
        neutral_zone = self._bfs(neutral_zone, band_non_neighbors)

        return neutral_zone

    def _bfs(self, neutral_zone, band_non_neighbors):
        """
        Returns the complete neutral zone by performing a breadth-first search
        of the neutral zone to check if any of the elements in
        score-band-non-neighbors are connected to any of the neutral zone
        members. The search is dynamic, i.e., as connected elements are found,
        these are themselves then considered members of the neutral zone, and
        used as seed elements for further search.

        :param neutral_zone: (list) of tuples of the form (genotype, score),
                which represents the neutral zone to use as a starting point for
                the BFS.
        :param band_non_neighbors: (list) of tuples of the form
                (genotype, score). This is the list of elements against which to
                check connectivity of elements of the given neutral zone.

        :return: (list) The extended neutral zone.

        """

        i = 0  # Index variable

        # Keep iterating for as long as there are elements left in the
        # neutral zone
        while i < len(neutral_zone):
            # Set of members of score-band-non-neighbors that are found to be
            # indirect neighbors of the focal element, and should therefore be
            # removed from the list of score-band-non-neighbors
            non_neighbors_to_remove = []

            # For each member of score-band-non-neighbors,
            for non_neighbor in band_non_neighbors:
                # Flag to indicate whether the non-neighbor is connected to the
                # ith element in the neutral zone,
                are_connected = self._are_neighbors(
                    neutral_zone[i]['sequence'],
                    non_neighbor['sequence']
                )

                # If the two are connected,
                if are_connected:
                    # Append the non-neighbor to the neutral zone
                    neutral_zone.append(non_neighbor)

                    # Mark the non-neighbor, so that it is removed from the
                    # set of score-band-non-neighbors after the for loop.
                    non_neighbors_to_remove.append(non_neighbor)

            # Remove those elements from the set of score-band-non-neighbors
            # that have already been added to the neutral zone
            band_non_neighbors = [
                e for e in band_non_neighbors
                if e not in non_neighbors_to_remove
            ]

            # Move on to the next element in the neutral zone
            i += 1

        return neutral_zone

    def _neighbors_processed_element(self, neutral_zone):
        """
        Checks if any member of the given neutral zone has a 1-neighbor outside
        the neutral zone that has already been processed (either as a member of
        a peak, or a member of a valley).

        :param neutral_zone: (list) of tuples of the form (genotype, score).

        :return: (bool) 'True' if any member of the neutral zone has an already
                    processed element as 1-neighbor . 'False' otherwise.

        """

        for element in neutral_zone:
            # All 1-neighbors of the element
            neighbors = self._get_neighbors_of(element['sequence'])

            # If any neighbor has already been processed,
            if not neighbors.isdisjoint(self._processed_genotypes):
                return True

        return False

    @staticmethod
    def _create_peak(neutral_zone):
        """
        Extracts all genotypes from the given neutral zone and returns them as
        a single collection that represents all members of a new peak.

        :param neutral_zone: (list) of tuples of the form (genotype, score).

        :return: (set) All genotypes in the peak.

        """

        return {e['sequence'] for e in neutral_zone}

    def _get_neighbors_of(self, genotype):
        """
        Returns all 1-neighbors of the given genotype.

        :param genotype: (str) Genotype for which to return the neighbors.

        :return: (set) All 1-neighbors of the given genotype.

        """

        return self._neighbor_map[genotype]

    def _are_neighbors(self, g_1, g_2):
        """
        Checks if g_2 is a 1-neighbor of g_1.

        :param g_1: (str) The first genotype.
        :param g_2: (str) The second genotype.

        :return: (bool) 'True' if g_2 is a 1-neighbor of g_1. 'False' otherwise.

        """

        return True if g_2 in self._get_neighbors_of(g_1) else False
