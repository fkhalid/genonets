"""
    peaks_functions
    ~~~~~~~~~~~~~~~

    Contains functions used for computation of peaks.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

from genonets_utils import Utils


# Provides functions to analyze the landscape defined by the geontype
# network, and detect peaks in the landscape.
# TODO: Once the peaks dict has been created for a network, perhaps the
#		path list can be deleted to release memeory ...
class PeakAnalyzer:
    # Constructor
    def __init__(self, network, netUtils, delta):
        # Reference to the network on which to perform this
        # analysis
        self.network = network

        # Get a reference to the NetworkUtils object
        self.netUtils = netUtils

        # Get a reference to the BitSeqManipulator in use
        self.bitManip = netUtils.bitManip

        # Keep a copy of the delta value
        self.delta = delta

        # Set once peaks have been identified
        self.peaks = None

    # Returns the peak dict which has seq as one of the keys
    def getPeakWithSeq(self, seq):
        # Get all peaks
        peaks = self.getPeaks()

        # For each peak
        for peakId in peaks.keys():
            # Get the peak dict
            peak = peaks[peakId]

            # Check if the given seq is in the peak
            if seq in peak["sequences"]:
                return peak

        # TODO: replace with proper exception handling ...
        print("Error: Count not find " + seq + " in any peak!!!!")
        return None

    # The function to be exposed to the calling object. Triggers the
    # peak detection algorithm and returns the result as a dictionary
    # with peak Ids as keys.
    def getPeaks(self, recompute=False):
        # If peak computation has not been done, or the caller has
        # explicitly asked for re-running the algorithms,
        if not self.peaks or recompute:
            # Get a structured array with tuples (sequence, escore),
            # sorted in descending order of e-scores.
            sortedArr = Utils.getSortedSeqEscArr(self.network,
                                                 self.bitManip.seqLength,
                                                 sortOrder="descending")

            # Identify peaks
            self.peaks = self.buildPeaks(sortedArr)

        return self.peaks

    # Orchistration function for the peak detection algorithm.
    def buildPeaks(self, elements):
        peaks = {}

        # List of elements that have been processed independently, or as
        # part of a neighborhood, e.g., all elements in a plateau.
        processed = []

        # Process each element, where elements are sorted in the
        # descending order of e-scores.
        for i in range(len(elements)):
            # If the element is in the list of already processed elements,
            if elements[i] in processed:
                # Skip iteration
                continue

            # Build a list of elements that belong to the neutral
            # zone, i.e., build the bin/band.
            binMembers = self.getBin(i, elements)

            # Build a list of elements within the bin that are
            # connected. This requires a BFS. If the iteration element
            # turns out to be a peak, all elements in this list will
            # be part of the plateau. In addition, get a list of elements
            # that are in the bin, but not connected.
            neutralZone, nonNeighs = self.getNeutralZone(elements[i], binMembers)
            # Append the element itself as well
            neutralZone.append(elements[i])

            # Get a list of peaks for which the focal element falls within
            # the bin.
            peakNonNeighs = self.getPeaksNonNeighs(elements[i], peaks)

            # Check connectivity with existing peaks. If the focal element
            # is connected to a peak, and it does not lie within the bin for
            # that peak, the entire neutral zone is a 'non-peak'.
            neighboringPeaks = self.getNeighPeaks(neutralZone, peakNonNeighs, peaks)

            # Depending on the results so far, append the list to an
            # existing peak, or create a new peak.
            if neighboringPeaks:
                # Append the entire neutral zone to the path list for the
                # neighboring peaks.
                self.appendToPeaks(neutralZone, neighboringPeaks, peaks)

                # Append the neutral zone to the processed list. If an entire
                # zone has been detected as connected to a peak, none of the
                # sequences can be peaks, and therefore do not need to
                # be processed independently.
                if len(neutralZone) > 1:
                    processed.extend(neutralZone)
            else:
                # Create a new peak, and add all members of the neutral zone
                # to this peak.
                peakId = len(peaks)
                peaks[peakId] = self.createPeak(neutralZone, nonNeighs)

                # If the peak is a plateau, mark all elements as processed,
                # so that none of them is processed independently.
                if len(neutralZone) > 1:
                    processed.extend(neutralZone)

        return peaks

    # Returns a list of elements that fall within the noise determined score
    # zone for the given focal element, identified by the given index.
    def getBin(self, index, elements):
        elBin = []

        focalElement = elements[index]

        # Starting from the next element in the list and until the
        # end
        for i in range(index + 1, len(elements)):
            # If the escore is within the band
            if elements[i]['escore'] >= focalElement['escore'] - self.delta:
                elBin.append(elements[i])
            else:
                # Once we've reached an element which lies in a lower
                # band, we know that the remaining elements will also
                # lie in lower bands. Therefore, we don't need to look
                # any further.
                break

        return elBin

    # The neutral zone consists of all elements that satisfy the following
    # conditions:
    # 1) The element lies within the bin for the given focal element
    # 2) The element is either a 1-neighbor of the focal element, or
    # 3) The element is an 'k'-neighbor of the focal element, such that all
    #	 elements that constitute the path between the element and the focal
    #	 element lie within the bin.
    # Returns the neutral zone, as well as a list of elements that are within
    # the bin, but are not connected to the focal element.
    def getNeutralZone(self, focalElement, binMembers):
        # Get the list of sequences that constitute the neighborhood of the
        # focal element
        neighSeqs = self.getNeighSeqs(focalElement)

        # Get a list of elements that exist in the neighborhood of the
        # focal element. These are by definition already part of the neutral
        # zone.
        neutralZone = [item for item in binMembers if item['sequence'] in neighSeqs]

        # Get all non-neighbor elements in the bin
        nonNeighs = [binMembers[i] for i in range(len(binMembers)) \
                     if binMembers[i] not in neutralZone]

        # Do a breadth-first search of the neutral zone to find elements in
        # the non-neighbors list that might be indirected connected.
        neutralZone, nonNeighs = self.bfsNeutralZone(neutralZone, nonNeighs)

        return neutralZone, nonNeighs

    # Performs a breadth-first search of the neutral zone to check if any of
    # the elements in the list of non-neighboring bin memebers are connected
    # to any of the members in the neutral zone. The search is dynamic, i.e.,
    # as connected elements are found, these are themselves then considered
    # as members of the neutral zone, and used as seed elements for further
    # search.
    # Returns the complete neutral zone, as well as a list of non-neighbors
    # left after the BFS.
    def bfsNeutralZone(self, neutralZone, nonNeighs):
        i = 0  # Index variable

        # Keep iterating for as long as there are elements left in the
        # neutral zone.
        while i < len(neutralZone):
            removeList = []

            # For each non-neighbor
            for nonNeigh in nonNeighs:
                # If the non-neighbor is connected to the ith element in the
                # neutral zone,
                if self.netUtils.areConnected(nonNeigh['sequence'], \
                                              neutralZone[i]['sequence']):
                    # Append the non-neighbor to the neutral zone
                    neutralZone.append(nonNeigh)

                    # Mark the non-neighbor, so that it is removed from the
                    # non-neighbors list after the for loop.
                    removeList.append(nonNeigh)

            # Remove the elements from the non-neighbors list that have already
            # been added to the neutral zone.
            nonNeighs = [item for item in nonNeighs if item not in removeList]

            # Move on to the next element in the neutral zone
            i += 1

        return neutralZone, nonNeighs

    # Checks if the given element was identified as a non-neighbor in the
    # bin for an existing peak.
    def getPeaksNonNeighs(self, focalElement, peaks):
        # Note: We do not check whether or not the e-score of the focal
        #		element + delta >= e-score of the peak summit, because
        #		if it weren't, the focal element here would not have
        #		made it to the non-neighbors list.
        return [peakId for peakId in peaks.keys() \
                if focalElement in peaks[peakId]["non-neighbors"]]

    # For each element in the neutral zone, checks whether it is connected
    # to an existing peak, except for peaks that are in 'peakNonNeighs'.
    # Returns a list of 'peakIds' for the neighboring peaks if any.
    def getNeighPeaks(self, neutralZone, peakNonNeighs, peaks):
        # List of peak Ids of all peaks to which this neutral zone is
        # connected.
        neighboringPeaks = []

        # For each element in the neutral zone
        for element in neutralZone:
            # Get peaks Id of all neighboring peaks for this element
            peaksAsNeighs = self.getNeighPeaksFor(element, peakNonNeighs, peaks)

            # The list should be empty if no existing peaks neighbor this
            # element
            if peaksAsNeighs:
                neighboringPeaks.extend(peaksAsNeighs)

        # Remove repeated values
        neighboringPeaks = list(set(neighboringPeaks))

        return neighboringPeaks

    # Checks if the given element is connected to a peak that is not in
    # 'peakNonNeighs'.
    # Returns a list of 'peakIds' for the neighboring peaks if any.
    def getNeighPeaksFor(self, element, peakNonNeighs, peaks):
        peakIds = []

        # For each existing peak,
        for peakId in peaks.keys():
            # If the focal element lies within the bin for this peak,
            if peakId in peakNonNeighs:
                # No need to check for indirect connectivity. Move on
                # to the next peak.
                continue

            # Get the list of connected elelments for this peak
            pathList = peaks[peakId]["list"]

            # For each element in the path list,
            for item in pathList:
                # Check if the given element is connected to current item
                # in the path list. The assumption here is that the first
                # items in the list are the peak sequences themselves.
                if self.netUtils.areConnected(item['sequence'], element['sequence']):
                    # The neutral zone has a connected to the peak. Therefore,
                    # add the peak ID to the list of connected peaks.
                    peakIds.append(peakId)

                    # Move on to the next peak
                    break

        return peakIds

    # Creates a peak which comprises the given neutral zone. Also, stores
    # the list of unconncted bin members.
    def createPeak(self, neutralZone, nonNeighs):
        # Dictionary: There are two keys, 1) List of sequences that
        # constitute the peak, and 2) list of sequences that are
        # connected to the peak.
        peak = {}

        # List sequences that constitute the peak. Reverse the neutral zone
        # so that the focal element that was appended at the end becomes
        # the head.
        neutralZone.reverse()
        peak["sequences"] = [neutralZone[i][0] \
                             for i in range(len(neutralZone))]

        # Initialize the list of the connected sequences with sequences
        # in the peak. This simplifies comparisons when searching through
        # the list for connectivity.
        peak["list"] = neutralZone

        # Store the non-neighboring elements from the same band to simplify
        # lookup operation during peak determination.
        peak["non-neighbors"] = nonNeighs

        return peak

    # Appends the given neutral zone to all the given neighboring peaks.
    def appendToPeaks(self, neutralZone, neighboringPeaks, peaks):
        # For each neighboring peak
        for peakId in neighboringPeaks:
            # Append the neutral zone to the list of elements found to
            # be connected to the peak.
            peaks[peakId]["list"].extend(neutralZone)

    # Returns the list of sequences that constitute the neighborhood of the
    # given element
    def getNeighSeqs(self, element):
        return [self.network.vs[neighbor]["sequences"] \
                for neighbor in self.netUtils.getNeighbors( \
                element['sequence'], self.network)]
