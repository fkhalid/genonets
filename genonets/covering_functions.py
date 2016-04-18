"""
    covering_functions
    ~~~~~~~~~~~~~~~~~~

    Contains functions used for computation of 'phenotype space covering'.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""


class CoveringAnalyzer:
    def __init__(self, network, rep_to_giant):
        # Reference to the network on which to perform this analysis
        self.network = network

        # Reference to dict: key=repertoire, value=giant
        self.repToGiant = rep_to_giant

    def covering_all(self):
        pass

    def covering(self, genotype):
        pass
