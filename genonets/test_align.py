#!/usr/bin/env python

from genonets_interface import Genonets
from cmdl_handler import CmdParser


def process(args):
    # Create the Genonets object. This will load the input file into
    # memory.
    gn = Genonets(args)

    # Use 'gn' to create genotype networks for all genotype sets.
    gn.create()


if __name__ == "__main__":
    # Parse the command line arguments using the Genonets command line handler, and
    # pass the list of arguments to 'process()'.
    process(CmdParser().getArgs())

    # Print message to indicate processing is done.
    print("Done.\n")
