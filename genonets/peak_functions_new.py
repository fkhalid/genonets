
import enum
import collections

from tqdm import tqdm

from genonets_utils import Utils


# GenotypeInfo = collections.namedtuple(
#     'GenotypeInfo',
#     ['neighbors', 'status']
# )
# GenotypeInfo.__doc__ = """
#     Represents an a record with important information about a genotype, e.g.,
#     a list of its neighbors, its current status, i.e., whether it has been
#     identified as a peak or non-peak, or if it has not been processed as at all.
# """


class PeakAnalyzer:
    # Constructor
    def __init__(self, network, net_utils, delta):
        # Reference to the network on which to perform this analysis
        self._network = network

        # Reference to the NetworkBuilder object
        self._net_utils = net_utils

        # Reference to the BitSeqManipulator in use
        self._bit_manip = net_utils.bitManip

        # Keep a copy of the delta value
        self._delta = delta

        # To hold results after computation
        self._peaks = None

        # Map {sequence: [all 1-neighbors]}
        self._neighbor_map = {
            s: set(self._net_utils.getNeighborSequences(s, self._network))
            for s in network.vs['sequences']
        }

        # Map {sequence: (category, i.e., peak or non-peak)}.
        # Keys are sequences that have been processed independently, or as
        # part of a neighborhood, e.g., all elements in a plateau. Values tell
        # use whether the sequences were identified as peaks or non-peaks.
        self._processed_sequences = {}

    def _get_neighbors_of(self, sequence):
        return self._neighbor_map[sequence]

    def _are_neighbors(self, sequence_1, sequence_2):
        return True if sequence_2 in self._neighbor_map[sequence_1] else False

    def _build_peaks(self, elements):
        """

        :param elements: (tuple) of the form (sequence, score).

        :return:

        """

        peaks = {}

        # Process each element, assuming the elements are sorted in the
        # descending order of scores
        for i in tqdm(xrange(len(elements))):
            # If the current element has already been processed,
            if elements[i]['sequence'] in self._processed_sequences:
                # Skip iteration
                continue

            # Build a list of elements for which the score lies within band
            # defined by the score of the current focal element
            score_band = self._build_score_band(i, elements)

            # Build a list of elements within the score band that are connected.
            # This requires a breadth-first search. If the focal element
            # turns out to be a peak, all elements in this list will be part of
            # the plateau.
            # In addition, get a list of elements that are in the score band,
            # but not connected to the neutral zone.
            neutral_zone, band_non_neighbors = \
                self._build_neutral_zone(elements[i], score_band)
            # Append the element itself as well
            neutral_zone.append(elements[i])

            # Get a list of peaks for which the focal element falls within
            # the bin.
            peakNonNeighs = self.getPeaksNonNeighs(elements[i], peaks)

            # Check connectivity with existing peaks. If the focal element
            # is connected to a peak, and it does not lie within the bin for
            # that peak, the entire neutral zone is a 'non-peak'.
            neighboringPeaks = self.getNeighPeaks(neutral_zone, peakNonNeighs, peaks)

            # Depending on the results so far, append the list to an
            # existing peak, or create a new peak.
            if neighboringPeaks:
                # Append the entire neutral zone to the path list for the
                # neighboring peaks.
                self.appendToPeaks(neutral_zone, neighboringPeaks, peaks)

                # Append the neutral zone to the processed list. If an entire
                # zone has been detected as connected to a peak, none of the
                # sequences can be peaks, and therefore do not need to
                # be processed independently.
                if len(neutral_zone) > 1:
                    self._processed_sequences.update({
                        e['sequence'] for e in neutral_zone
                    })
            else:
                # Create a new peak, and add all members of the neutral zone
                # to this peak.
                peakId = len(peaks)
                peaks[peakId] = self.createPeak(neutral_zone, band_non_neighbors)

                # If the peak is a plateau, mark all elements as processed,
                # so that none of them is processed independently.
                if len(neutral_zone) > 1:
                    self._processed_sequences.update({
                        e['sequence'] for e in neutral_zone
                    })

        return peaks

    def _build_score_band(self, index_of_focal, elements):
        """
        Returns a list of elements that fall within the noise determined score 
        band for the given focal element, where the focal element is identified
        by the given index.

        :param index_of_focal: (int) Index of the focal element in the given
                list of all elements.
        :param elements: (numpy.array) Array of (sequence, score) tuples
                corresponding to all genotypes in the network in consideration.
        
        :return: (list) All elements for which the score lies within the band 
                 of the focal element.
        
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
                if elements[i]['sequence'] in self._processed_sequences:
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
        The neutral zone consists of all elements that satisfy the following
        conditions:
            1) The element lies within the score band for the given focal 
               element.
            2) The element is either a 1-neighbor of the focal element, or
            3) The element is a 'k'-neighbor of the focal element, such that all
               elements that constitute the path between the element and the 
               focal element lie within the score band of the focal element.
        
        Returns the neutral zone, as well as a list of elements that are within
        the score band, but are not connected to the focal element.
        
        :param focal_element: (tuple) of form (sequence, score) corresponding
                to the element for which the score band was constructed.
        :param score_band: (list) of (sequence, score) tuples, where each
                tuple represents a sequence for which the score lies in the
                score band for the focal element.
        
        :return: (tuple) (
                            (list) neutral zone,
                            (list) element in the given score band that are
                            not 1-neighbors of the focal element
                        )
        
        """

        # Sequences that constitute the neighborhood of the focal element
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
        band_non_neighbors = set(score_band) - set(neutral_zone)

        # Do a breadth-first search of the neutral zone to find those elements
        # in band_non_neighbors that might be indirectly connected (via other
        # members of the neutral zone) to the focal element.
        neutral_zone, band_non_neighbors = \
            self._bfs(neutral_zone, band_non_neighbors)

        return neutral_zone, band_non_neighbors

    def _bfs(self, neutral_zone, band_non_neighbors):
        """
        Performs a breadth-first search of the neutral zone to check if any of
        the elements in score-band-non-neighbors are connected to any of the
        neutral zone members. The search is dynamic, i.e., as connected
        elements are found, these are themselves then considered members of the
        neutral zone, and used as seed elements for further search.

        Returns the complete neutral zone, as well as score-band-non-neighbors
        left after the BFS.

        :param neutral_zone:
        :param band_non_neighbors:

        :return:

        """

        i = 0  # Index variable

        # Keep iterating for as long as there are elements left in the
        # neutral zone.
        while i < len(neutral_zone):
            # Set of members of score-band-non-neighbors that are found to be
            # indirect neighbors of the focal element, and should therefore be
            # removed from the list of score-band-non-neighbors.
            non_neighbors_to_remove = set()

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
                    non_neighbors_to_remove.add(non_neighbor)

            # Remove those elements from the set of score-band-non-neighbors
            # that have already been added to the neutral zone
            band_non_neighbors -= non_neighbors_to_remove

            # Move on to the next element in the neutral zone
            i += 1

        return neutral_zone, band_non_neighbors

    # Checks if the given element was identified as a non-neighbor in the
    # bin for an existing peak.
    def getPeaksNonNeighs(self, focal_element, peaks):
        # Note: We do not check whether or not the e-score of the focal
        #		element + delta >= e-score of the peak summit, because
        #		if it weren't, the focal element here would not have
        #		made it to the non-neighbors list.
        return [
            peak_id for peak_id in peaks
            if focal_element['sequence'] in peaks[peak_id]['non-neighbors']
        ]
