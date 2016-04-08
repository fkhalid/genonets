"""
    graph_utils
    ~~~~~~~~~~~

    Function for network manipulation.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import json
import igraph

from genonets.seq_align import Aligner

# Provides methods for constructing and manipulating genotype
# networks. The methods in this class are agnostic to molecule
# type.
class NetworkBuilder:
    # Constructor
    def __init__(self, seqBitManip):
        # Store reference to the sequence/bit manipulator object
        # in use
        self.bitManip = seqBitManip

        # TODO: Initialize the aligner object here ...
        # FIXME: Figure out a way to access sequence length and alphabet ...
        self.aligner = Aligner(alphabet="DNA",
                               seqlen=8,
                               match=1,
                               mismatch=-1,
                               gap=-0.5)

    # Return 'True' if the source sequence is connected to the target
    # sequence, and 'False' otherwise.
    def areConnected(self, source, target):
        return self.bitManip.areNeighbors(self.bitManip.seqToBits(source),
                                          self.bitManip.seqToBits(target))

    # Create a genotype network using the give list of sequences
    def createGenoNet(self, repertoire, sequences, scores):
        # Initialize the empty network
        network = igraph.Graph()

        # Set repertoire name as a graph attribute
        network["name"] = repertoire

        # Create vertices
        network.add_vertices(len(sequences))

        # Add scores to all vertices
        # FIXME: Change back from 'escores' to 'scores'
        network.vs["escores"] = scores

        # Add sequences to all vertices
        network.vs["sequences"] = sequences

        # # Add bit formatted sequences to all vertices
        # # Note: For very long sequences, the size of the corresponding
        # # integer value becomes too long for igraph to be able
        # # to deal with it (since igraph is C based). Therefore,
        # # bitseqs cannot be stored in the graph.
        # bitseqs = [self.bitManip.seqToBits(seq) for seq in sequences]
        #
        # # Get a list of pairs of indices which represent paris of
        # # sequences that are 1-neighbors, i.e., only one mutation
        # # apart.
        # edges = [
        #     (i, j)
        #     for i in range(len(bitseqs) - 1)
        #     for j in range(i + 1, len(bitseqs))
        #     if self.bitManip.areNeighbors(bitseqs[i], bitseqs[j])
        # ]

        # Get a list of pairs of indices which represent paris of
        # sequences that are 1-neighbors, i.e., only one mutation
        # apart.
        edges = [
            (i, j)
            for i in range(len(sequences) - 1)
            for j in range(i + 1, len(sequences))
            if self.aligner.areNeighbors(sequences[i],
                                         sequences[j],
                                         self.bitManip)
        ]

        # Connect the two with an edge
        network.add_edges(edges)

        return network

    # For the given sequence, returns a list of 1-neighbors that are not part of
    # the given genotype network.
    # Note: The list of external neighbors returned comprises bitseqs, not
    # strings. The caller must be aware of this.
    def getExternalNeighbors(self, sequence, network):
        # Get a list of all possible 1-neighbors of the given sequence
        allNeighbors = [
            neighbor for neighbor in
            self.bitManip.generateNeighbors(self.bitManip.seqToBits(sequence))
        ]

        # Get a list of 1-neighbors that are part of the given genotype
        # network
        genoNetNeighbors = [
            self.bitManip.seqToBits(network.vs[neighbor]["sequences"])
            for neighbor in self.getNeighbors(sequence, network)
        ]

        # Get a list of sequences that are 1-neighbors, but not part of the
        # genotype network. This is achieved by taking the difference of
        # two sets, i.e., elements in allNeighbors not in genoNetNeighbors
        externNeighbors = list(set(allNeighbors) - set(genoNetNeighbors))

        return externNeighbors

    # Get all unique external neighbors for the given genotype network
    def getAllExtNeighbors(self, network):
        # List to be populated with external neighbors
        extNeighbors = []

        # Get all sequences
        sequences = network.vs["sequences"]

        # For each sequence, get the list of external neighbors and
        # append it to extNeighbors
        for sequence in sequences:
            extNeighbors.extend(self.getExternalNeighbors(sequence, network))

        # Remove redundant sequences
        extNeighbors = list(set(extNeighbors))

        return extNeighbors

    # Get the vertex based on the sequence string from the
    # given network
    def getVertex(self, sequence, network):
        try:
            return network.vs.find(sequences=sequence)
        except ValueError:
            print("Error! ... Non-existent vertex requested: " + str(sequence))

    # Returns a list of sequences that constitute the neighborhood
    # of the given sequence within the given network.
    def getNeighbors(self, sequence, network):
        return network.neighbors(self.getVertex(sequence, network))

    # Plot the given network using the given layout.
    def plotNetwork(self, network, layout, outPath):
        vertices = [network.vs[i].index for i in range(0, network.vcount())]
        # For GNs,
        # degrees = network.degree(vertices)
        # For ENs,
        degrees = network.outdegree(vertices)

        visual_style = {}
        # For GNs:
        # visual_style["vertex_size"] = [x*2 for x in degrees]
        # For ENs:
        # visual_style["vertex_size"] = [x*2 for x in degrees]
        visual_style["vertex_size"] = [
            (x * 10 * x * 10)
            for x in network.vs["Evolvability"]
        ]
        # visual_style["edge_width"] = [math.log(x) for x in network.es["weight"]]
        # visual_style["vertex_color"] = [int(x*300) for x in network.vs["escores"]]

        layout = network.layout(layout)

        igraph.plot(network, outPath + network["name"] + "_evonet.svg", bbox=(1500, 1500), **visual_style)

    # Return a list, the size of which represents the number of commnected
    # component in the network. Each element in the list corresponds to the
    # size of one of the connected components.
    def getComponents(self, network):
        vertexCluster = network.components()

        return vertexCluster.sizes()

    # Returns the gian component.
    def getGiantComponent(self, network):
        return network.components().giant()

    # Create a network that consists of the given sequence and its
    # one-neighbors.
    def createNeighborhoodNet(self, sequence, network):
        sequences = [sequence]
        escores = []

        # Get the vertex corresponding to the sequence
        vertex = self.getVertex(sequence, network)

        # Get all adjacent vertices
        neighbors = network.neighbors(vertex)

        # Get escore for the given sequence
        escores.append(network.vs[vertex.index]["escores"])

        # Get sequences for all neighbors
        for neighbor in neighbors:
            sequences.append(network.vs[neighbor]["sequences"])
            escores.append(network.vs[neighbor]["escores"])

        # Initialize the empty network
        neighborhoodNet = igraph.Graph()

        # Create the neighborhood network

        # Create vertices
        neighborhoodNet.add_vertices(len(sequences))

        # Add attributes
        neighborhoodNet.vs["label"] = [
            sequences[i] + "\n" + str(escores[i])
            for i in range(len(sequences))
        ]

        for i in range(1, len(sequences)):
            neighborhoodNet.add_edge(0, i)

        # Print the list of sequences for validation purposes
        print(sequences)

        return neighborhoodNet

    # Create a network of the given repertoires where repertoire 'A' is
    # connected to repertoire 'B' if 'A' has one or more external
    # neighbors in 'B'.
    # Arguments: title=network name, genoNets=giant for the repertoire
    # TODO: We have the option of either using the GNs as nodes, or just
    # the GN names. I need to think more about the pros and cons
    # of each approach ...
    def createEvoNet(self, title, genoNets):
        # Get a list of repertoire names
        repertoires = [
            gn["name"].rpartition("_dominant")[0]
            if gn["name"].rpartition("_dominant")[0]
            else gn["name"]
            for gn in genoNets
        ]
        # repertoires = [gn["name"].replace("_dominant", "") for gn in genoNets]
        # repertoires = [gn["name"].rstrip("_dominant") for gn in genoNets]

        # Get a list of evolvability scores corresponding to the
        # repertoires
        evolvability = [gn["Evolvability"] for gn in genoNets]

        # Initialize the empty network
        network = igraph.Graph(directed=True)

        # Set title as the graph name
        network["name"] = title

        # Create vertices: Each repertoire name is a vertex in the network
        network.add_vertices(len(repertoires))

        # Add attributes

        # TODO: Remove this from the tool ...
        labels = [
            repertoires[i] + "\n" + str(float("{0:.2f}".format(evolvability[i])))
            for i in range(len(repertoires))
        ]
        # labels = []
        # labels[:] = repertoires

        # Add labels to all vertices
        network.vs["label"] = labels  # TODO: Remove this from the tool ...
        # Add repertoire names to all vertices
        network.vs["Repertoires"] = repertoires
        network.vs["GenotypeSet"] = repertoires  # Redundant; added for GML
        network.vs["Evolvability"] = evolvability

        # Get a list of pairs of indices which represent paris of
        # repertoires that are neighbors, i.e., share at least one external
        # neighbor. The edge determination algorithm takes into account the
        # antisymmetric relation.
        edges = [
            (i, j)
            for i in range(len(repertoires))
            for j in range(len(repertoires))
            if i != j and self.repertoiresAreConnected(genoNets[i], repertoires[j])
        ]

        # Add the above computed edges to the network
        network.add_edges(edges)

        # Add edge weights

        # Create a dict: keys=repertoires, values=genonets
        repToNetDict = {
            repertoires[i]: genoNets[i] for i in range(len(repertoires))
        }

        # Set weight for each edge, where weight=No. of evo target genotypes from i to j
        weights = [
            self.getEdgeWeight(edges[i], network, repToNetDict)
            for i in range(len(edges))
        ]

        # Add edge weights and edge labels as attributes
        network.es["weight"] = weights
        network.es["label"] = weights

        return network

    # Returns True if 'gn1' has at least one external neighbor in the GN for
    # 'rep2'
    def repertoiresAreConnected(self, gn1, rep2):
        # Result flag
        areConnected = False

        # Get the list of evolvability target repertoires for gn1
        evoTargets = json.loads(gn1["Evolvability_targets"])

        # If there is at least one target in the list,
        if evoTargets:
            # If 'rep2' is in the list of targets,
            if rep2 in evoTargets:
                # 'gn1' has external neighbors in rep2
                areConnected = True

        return areConnected

    # Return the weight for the given edge of the given EN
    def getEdgeWeight(self, edge, en, repToNetDict):
        # Get the repertoire names for source and target vertices
        srcRep = en.vs[edge[0]]["Repertoires"]
        trgRep = en.vs[edge[1]]["Repertoires"]

        # Get the source GN
        srcGN = repToNetDict[srcRep]

        # Get the list of target dicts in source GN
        allTargets = srcGN.vs["Evolvability_targets"]

        # External neighbors of source in target; to be
        # populated in the loop below.
        extNeighbors = []

        # Extract the required dict
        for trgDict in allTargets:
            # If the target repertoire is a key in this dict,
            if trgRep in trgDict:
                # Add external neighbors to the global list
                extNeighbors.extend(trgDict[trgRep])

        # Remove redundant genotypes
        extNeighbors = list(set(extNeighbors))

        return len(extNeighbors)
