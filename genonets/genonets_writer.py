"""
    genonets_writer
    ~~~~~~~~~~~~~~~

    Helper functions for writing results to files.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import os
import copy
import warnings

from genonets_constants import ErrorCodes
from genonets_exceptions import GenonetsError
from genonets_constants import GenonetsConstants as gc


class Writer:
    @staticmethod
    def writeInParamsToFile(paramsDict, path):
        fileName = path + "in_params.txt"

        # If the required directories have not already been created
        if not os.path.exists(os.path.dirname(path)):
            # Create the directories
            try:
                os.makedirs(os.path.dirname(path))
            except os.error:
                print("Error: " +
                      ErrorCodes.getErrDescription(ErrorCodes.CANNOT_CREATE_DIRECTORY) +
                      ": Path - " + path)

                raise GenonetsError(ErrorCodes.CANNOT_CREATE_DIRECTORY,
                                    "Path - " + path)

        # Open the file
        try:
            with open(fileName, "w") as outFile:
                # For each input parameter,
                for param in paramsDict.keys():
                    outFile.write(param + ": " + paramsDict[param] + "\n")
        except Exception:
            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.CANNOT_WRITE_TO_FILE) +
                  ": Path - " + path)

            raise GenonetsError(ErrorCodes.CANNOT_WRITE_TO_FILE)

    @staticmethod
    def writeNetsToFile(repToNetDict, repToGiantDict, netBuilder, path, attrsToIgnore, repertoires=gc.ALL):
        # If all repertoires should be considered,
        if repertoires == gc.ALL:
            # Get a list of all repertoires
            repertoires = repToNetDict.keys()

        # If the required directories have not already been created
        if not os.path.exists(os.path.dirname(path)):
            # Create the directories
            os.makedirs(os.path.dirname(path))

        # For each repertoire
        for repertoire in repertoires:
            # Write the entire network to file in GML format
            Writer.writeNetToFile(repToNetDict[repertoire], path, attrsToIgnore)

            # Get the number of components in the network
            numComponents = len(netBuilder.getComponents(repToNetDict[repertoire]))

            # If the network as more than one components,
            if numComponents > 1:
                #  Write the giant component to the file in GML format
                Writer.writeNetToFile(repToGiantDict[repertoire], path, attrsToIgnore)

    @staticmethod
    def writeNetToFile(network, path, attrsToIgnore):
        # File name
        fileName = path + network["name"] + ".gml"

        # Create a deep copy of the network
        netToWrite = copy.deepcopy(network)

        # Get a list of network attributes
        netAttrs = netToWrite.attributes()

        # Remove the network attributes that need not be written to file
        for attr in attrsToIgnore(level="network"):
            if attr in netAttrs:
                del netToWrite[attr]

        # Get a list of vertex attributes
        vtxAttrs = netToWrite.vs.attributes()

        # Remove the vertex attributes that need not be written to file
        for attr in attrsToIgnore(level="vertex"):
            if attr in vtxAttrs:
                del netToWrite.vs[attr]

        # If the given network is a genotype network, and not a
        # phenotype network,
        if not network.is_directed():
            # Rename 'sequences' and 'escores'
            netToWrite.vs["genotype"] = netToWrite.vs["sequences"]
            netToWrite.vs["score"] = netToWrite.vs["escores"]
            del netToWrite.vs["sequences"]
            del netToWrite.vs["escores"]

        # Write network to file in GML format, while ignoring 'igraph'
        # related warnings.
        with warnings.catch_warnings():
            # Create a filter for 'igraph' warnings
            warnings.simplefilter("ignore")

            # Write GML file
            netToWrite.write(fileName, format="gml")

    @staticmethod
    def writeNetAttribs(repToNetDict, repToGiantDict, netBuilder, path, attrsToIgnore, order, repertoires=gc.ALL):
        # If all repertoires should be considered,
        if repertoires == gc.ALL:
            # Get a list of all repertoires
            repertoires = repToNetDict.keys()

        # Construct the list of attributes to be written to file
        attributes = repToNetDict[repertoires[0]].attributes()
        attributes.remove("name")

        # Get the number of components in the network
        numComponents = len(netBuilder.getComponents(repToNetDict[repertoires[0]]))

        # If the network as more than one components,
        if numComponents > 1:
            # Get attributes from giant as well
            attributes.extend(repToGiantDict[repertoires[0]].attributes())

        # Remove the attributes that need not be written to file
        for attr in attrsToIgnore():
            if attr in attributes:
                attributes.remove(attr)

        # Sort the attributes, so that they are written to file in a
        # fixed, pre-determined order
        attributes.sort(key=order)

        # If the required directories have not already been created
        if not os.path.exists(os.path.dirname(path)):
            # Create the directories
            os.makedirs(os.path.dirname(path))

        fileName = path + "Genotype_set_measures.txt"
        dataFile = open(fileName, 'w')

        # Write the row of headers
        dataFile.write("Genotype_set" + "\t")
        for attribute in attributes:
            dataFile.write(attribute + "\t")

        dataFile.write("\n")

        # For each repertoire,
        for repertoire in repertoires:
            # Write the network name
            dataFile.write(repToNetDict[repertoire]["name"] + "\t")
            # For each attribute,
            for attribute in attributes:
                # If the attribute is in the main network,
                if attribute in repToNetDict[repertoire].attributes():
                    dataFile.write(str(repToNetDict[repertoire][attribute]) + "\t")
                else:
                    dataFile.write(str(repToGiantDict[repertoire][attribute]) + "\t")

            dataFile.write("\n")

        dataFile.close()

    @staticmethod
    def writeSeqAttribs(repToGiantDict, path, attrsToIgnore, order, repertoires=gc.ALL):
        # If all repertoires should be considered,
        if repertoires == gc.ALL:
            # Get a list of all repertoires
            repertoires = repToGiantDict.keys()

        # If the required directories have not already been created
        if not os.path.exists(os.path.dirname(path)):
            # Create the directories
            os.makedirs(os.path.dirname(path))

        for repertoire in repertoires:
            # Create the file name
            fileName = path + repertoire + "_genotype_measures.txt"

            # Open file to write
            with open(fileName, 'w') as dataFile:
                Writer.writeSeqAttribsFor(repToGiantDict[repertoire],
                                          dataFile, attrsToIgnore, order)

    @staticmethod
    def writeSeqAttribsFor(network, dataFile, attrsToIgnore, order):
        # Get vertex level attributes
        attributes = network.vs.attributes()

        # Remove the attributes that need not be written to file
        for attr in attrsToIgnore():
            if attr in attributes:
                attributes.remove(attr)

        # Sort the attributes, so that they are written to file in a
        # fixed, pre-determined order
        attributes.sort(key=order)

        # Write the row of headers
        dataFile.write("Sequence" + "\t")
        for attribute in attributes:
            dataFile.write(attribute + "\t")

        dataFile.write("\n")

        sequences = network.vs["sequences"]

        # For each sequence,
        for i in range(len(sequences)):
            # Write the sequence
            dataFile.write(sequences[i] + "\t")
            # For each attribute,
            for attribute in attributes:
                # Write value to file
                dataFile.write(str(network.vs[i][attribute]) + "\t")

            dataFile.write("\n")

    @staticmethod
    def writeOverlapToFile(overlapMat, repertoires, path):
        fileName = path + "Genotype_set_overlap.txt"

        # If the required directories have not already been created
        if not os.path.exists(os.path.dirname(path)):
            # Create the directories
            os.makedirs(os.path.dirname(path))

        dataFile = open(fileName, 'w')

        # Write header: All repertoire names as column headers
        for repertoire in repertoires:
            dataFile.write("\t" + repertoire)

        # Write the matrix with each row beginning with the
        # corresponding repertoire name
        for i in range(len(overlapMat)):
            dataFile.write("\n")

            # Write RBP name
            dataFile.write(repertoires[i])

            for j in range(len(overlapMat)):
                if not i == j:
                    dataFile.write("\t" + str(overlapMat[i][j]))
                else:
                    dataFile.write("\t" + "NaN")

        dataFile.close()
