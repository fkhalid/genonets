Genonets
========

This package provides a high level interface for construction and analysis of genotype networks from data. Also,
this is the Python package used by the `Genonets Server <http://ieu-genonets.uzh.ch/>`_.

Documentation, including tutorials and API documentation, is available `here <http://ieu-genonets.uzh.ch/python_package>`_.

----

New in version 1.1.5
~~~~~~~~~~~~~~~~~~~~

Bug fix: The fix affects the results of Robustness analysis; only when Genonets is used with '-rc' or
'--use_reverse_complements' options. The impact of this change is higher on genotype level results, but
minimal on genotype set level results. The details of the issue can be found
`here <https://github.com/fkhalid/genonets/issues/10>`_.

New in version 1.1.3
~~~~~~~~~~~~~~~~~~~~

Bug fix: The fix affects the results of Evolvability, Accessibility, Neighbor abundance, Diversity index, and Overlap analyses,
only when Genonets is used with '-rc' or '--use_reverse_complements' options. The impact of this change is higher on genotype
level results, but minimal on genotype set level results. The details of the issue can be found
`here <https://github.com/fkhalid/genonets/issues/9>`_.

New in version 1.1.0
~~~~~~~~~~~~~~~~~~~~

The public interface in ``genonets.genonets_interface.Genonets`` has been changed, i.e., several method signatures
used in the previous versions are no longer valid. Please see the API documentation `here <http://ieu-genonets.uzh.ch/python_package>`_.

New in version 1.0.7
~~~~~~~~~~~~~~~~~~~~

- An optional command line argument, '-v' or '--verbose' has been introduced. This enables the verbose mode. When used with python '-u' flag,  detailed progress information is printed to the standard output.
- A new analysis type 'PATHS_RATIOS' has been added. It enables the computation of ratio of 'accessible mutational paths' to 'all shortest mutational paths' for a given distance from summit.
- An optional command line argument, '-rc' or '--use_reverse_complements' has been introduced. This option can only be used with alphabet type 'DNA'. If this option is given, in addition to the genotypes, reverse complements of the genotypes are also considered during genotype network creation, as well as during 'Evolvability', 'Accessibility', 'Neighbor abundance', and 'Diversity index' analysis types.

Installation
------------

Linux (tested on Ubuntu 14.04 LTS and Ubuntu 16.04 LTS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Mac OS X El Capitan
~~~~~~~~~~~~~~~~~~~

We highly recommend using ``virtualenv``, or better yet, ``Anaconda``, for installation on Mac OS X El Capitan.

In case you do not already have ``virtualenv`` installed on your system, use the following command to install ``virtualenv``:

``pip install virtualenv``

In the directory of your choice, create a virtual environment. In the following example, we will create a virtual environment called ``venv_genonets``:

``virtualenv venv_genonnets``

Now, activate ``venv_genonets`` as follows:

``source venv_genonets/bin/activate``

You are now ready to install Genonets. Use the following command:

``pip install genonets``

Note: Every time you need to use ``genonets``, you will have to activate the corresponding virtual environment.

Windows
~~~~~~~

Instructions for Windows are basically the same, except in certain cases installation of dependencies fails. If that happens, follow these steps:

1. Download the 'whl' files for ``numpy`` and ``python-igraph`` from http://www.lfd.uci.edu/~gohlke/pythonlibs/. E.g.,

 i. ``numpy-1.10.2+mkl-cp27-none-win32.whl``
 ii. ``python_igraph-0.7.1.post6-cp27-none-win32.whl``

3. ``pip install python_igraph-0.7.1.post6-cp27-none-win32.whl``
4. ``pip install numpy-1.10.2+mkl-cp27-none-win32.whl``
5. And finally, ``pip install genonets``
