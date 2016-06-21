Genonets
========

This is the Python package used by the `Genonets Server
<http://ieu-genonets.uzh.ch/>`_ for creating and analyzing genotype networks from raw data. The details of the analyses used and the attributes computed can be found on the `Learn Genonets
<http://ieu-genonets.uzh.ch/learn>`_ page.

----

New in version 1.0.7
~~~~~~~~~~~~~~~~~~~~

- An optional command line argument, '-v' or '--verbose' has been introduced. This enables the verbose mode. When used with python '-u' flag,  detailed progress information is printed to the standard output.
- A new analysis type 'PATHS_RATIOS' has been added. It enables the computation of ratio of 'accessible mutational paths' to 'all shortest mutational paths' for a given distance from summit.
- An optional command line argument, '-rc' or '--use_reverse_complements' has been introduced. This option can only be used with alphabet type 'DNA'. If this option is given, in addition to the genotypes, reverse complements of the genotypes are also considered during genotype network creation, as well as during 'Evolvability', 'Accessibility', 'Neighbor abundance', and 'Diversity index' analysis types.

Installation
------------

Linux and Mac OS
~~~~~~~~~~~~~~~~

Using ``pip``,

``pip install genonets``

In case you get a 'permission' related error, try the following:

``sudo pip install genonets``

You can also install Genonets directly from the source package.

``python setup.py install``

Again, in case you run into permission related errors,

``sudo python setup.py install``

When trying to install genonets on a machine with ``Ubuntu 14.04 LTS`` that does not already have the required version of ``python-igraph`` installed, ``pip`` sometimes fails to install the C core of igraph. If that happens, follow these steps:

1. ``sudo apt-get install build-essential``
2. ``sudo apt-get python-dev``
3. ``sudo apt-get install libxml2-dev``
4. ``sudo apt-get install libz-dev``
5. ``sudo pip uninstall genonets``
6. Finally, ``sudo pip install genonets``

Windows
~~~~~~~

Instructions for Windows are basically the same, except in certain cases installation of dependencies fails. If that happens, follow these steps:

1. Download the 'whl' files for ``numpy`` and ``python-igraph`` from http://www.lfd.uci.edu/~gohlke/pythonlibs/. E.g.,
    i. ``numpy-1.10.2+mkl-cp27-none-win32.whl``
    ii. ``python_igraph-0.7.1.post6-cp27-none-win32.whl``

3. ``pip install python_igraph-0.7.1.post6-cp27-none-win32.whl``
4. ``pip install numpy-1.10.2+mkl-cp27-none-win32.whl``
5. And finally, ``pip install genonets``

Genonets quick start
--------------------

The best way to get started is to work through ``genonets_exmpl_simple.py`` available in the ``genonets/genonets/sample`` directory. In case you cannot locate the directory in which genonets is installed, you can download the source tarball from the genonets `PyPI page <https://pypi.python.org/pypi/genonets>`_, and find the ``sample`` folder inside the extracted ``genonets`` directory.

To get started, first copy the ``sample`` folder in a directory with write privileges. Then, try the following command:

``python genonets_exmpl_simple.py DNA true data/genonets_sample_input.txt 0.35 results_simple``

This command does the following:

1. Parses the sample input file located in the ``data`` directory
2. Creates genotype networks for all available genotype sets
3. Performs all available analyses on the genotype sets
4. Writes the following in the in the ``results_simple`` directory:
    i. A file with genotype network level attributes for all genotype sets
    ii. For each genotype network, a file with genotype level attributes
    iii. GML files for genotype networks

The following command can be used to view a description of command line arguments:

``python genonets_exmpl_simple.py -h``

The ``genonets/genonets/sample`` directory also includes other sample files, each highlighting different features.
