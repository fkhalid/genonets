
"""
    genonets_utils
    ~~~~~~~~~~~~~~

    General utility functions.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import numpy as np


class Utils :
    # Reverses the given dictionary, i.e., sets keys as values and
    # vice versa. Returns the reverse dictionary.
    # Note: This can only work if the values in the given dict are
    #		hashable. This function will not work for lists, etc.
    #		Keys and values must have a one-to-one mapping.
    @staticmethod
    def reverseDict(inDict) :
        return {inDict[key]: key for key in inDict.keys()}

    @staticmethod
    def getSeqWithMaxScore(network, seqLength) :
        # Get structured array with tuples: (sequence, escore),
        # sorted in ascending order of escores.
        sortedArr = Utils.getSortedSeqEscArr(network, seqLength,
                        sortOrder="ascending")

        # Last sequence in the array is the one with the highest
        # e-score.
        return sortedArr[-1]['sequence']

    # Returns a structured array with tuples: (sequence, escore),
    # sorted in either 'ascending' or 'descending' order of escores;
    # as determined by the 'sortOrder' received as argument.
    @staticmethod
    def getSortedSeqEscArr(network, seqLength, sortOrder) :
        # Get the structred array
        seqEscrArr = Utils.getSeqEscrArr(network, seqLength)

        # Sort the array in ascending order of escores
        sortedArr = np.sort(seqEscrArr, order="escore")

        # If the required order is 'ascending',
        if sortOrder == "ascending" :
            # We already have what we need
            return sortedArr

        # If the required order is 'descending',
        if sortOrder == "descending" :
            # Get a reverse view, which is equivalent of sorting
            # in descending order.
            revArr = sortedArr[::-1]

            return revArr

    # Generates a structured array with tuples: (sequence, escore) for
    # all genotypes in the network.
    @staticmethod
    def getSeqEscrArr(network, seqLength) :
        # Get a list of all sequences in the network
        seqs = network.vs["sequences"]

        # Get a list of all escores
        escores = network.vs["escores"]

        # Create the structured array
        dtype = [('sequence', 'S' + str(seqLength)),
                    ('escore', float)]
        values = [(seqs[i], escores[i]) for i in range(len(seqs))]

        return np.array(values, dtype=dtype)
