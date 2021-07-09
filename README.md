# Genonets

This package provides a high level interface for construction and analysis of genotype networks 
from data. It was also used as the analysis backend for the 
[Genonets Server](#Genonets-Server-end-of-life).

**Note:** Please [email us](mailto:genonets@outlook.com) for questions, comments, 
suggestions, and/or feature requests.

----

### New in version 1.1.10

New feature: If the `EPISTASIS` analysis is requested, the `Genotype_set_measures.txt` file now contains an additional 
column titled `Epistasis squares`. Data format for the new column and other detail are available 
[here](https://github.com/fkhalid/genonets/issues/18).

### New in version 1.1.9

Bug fix: In light of the issue reported [here](https://github.com/fkhalid/genonets/issues/16), the algorithm for
identification of peaks has been revised, as well as thoroughly validated and tested.


### New in version 1.1.8

Bug fix: In certain cases, the computation of peaks would result in two or more peaks that share
one or more genotypes. This behavior was incorrect, and has therefore been fixed. The corresponding
issue is reported [here](https://github.com/fkhalid/genonets/issues/14).

### New in version 1.1.7

Performance optimization: The code for computation of peaks has been optimized
so that it runs significantly faster than the previous version. The changes made
only affect performance; the algorithm remains the same. The corresponding issue is reported
[here](https://github.com/fkhalid/genonets/issues/12).

### New in version 1.1.6

Enhancement: The order in which genotype set names appear in the result file `Genotype_set_ovelrap.txt`, is now the
same order in which genotype set names appear in the input file. The corresponding issue is reported
[here](https://github.com/fkhalid/genonets/issues/11).

### New in version 1.1.5

Bug fix: The fix affects the results of Robustness analysis; only when Genonets is used with `-rc` or
`--use_reverse_complements` options. The impact of this change is higher on genotype level results, but
minimal on genotype set level results. The details of the issue can be found
[here](https://github.com/fkhalid/genonets/issues/10).

### New in version 1.1.3

Bug fix: The fix affects the results of Evolvability, Accessibility, Neighbor abundance, Diversity index, and Overlap 
analyses, only when Genonets is used with `-rc` or `--use_reverse_complements` options. The impact of this change is 
higher on genotype level results, but minimal on genotype set level results. The details of the issue can be found
[here](https://github.com/fkhalid/genonets/issues/9).

### New in version 1.1.0

The public interface in `genonets.genonets_interface.Genonets` has been changed, i.e., several method signatures
used in the previous versions are no longer valid. Please see the API documentation 
[here](http://ieu-genonets.uzh.ch/python_package).

### New in version 1.0.7

* An optional command line argument, `-v` or `--verbose` has been introduced. This enables the verbose mode. When used 
with python `-u` flag,  detailed progress information is printed to the standard output.
* A new analysis type `PATHS_RATIOS` has been added. It enables the computation of ratio of 
"accessible mutational paths" to "all shortest mutational paths" for a given distance from summit.
* An optional command line argument, `-rc` or `--use_reverse_complements` has been introduced. This option can only be 
used with alphabet type 'DNA'. If this option is given, in addition to the genotypes, reverse complements of the 
genotypes are also considered during genotype network creation, as well as during 'Evolvability', 'Accessibility', 
'Neighbor abundance', and 'Diversity index' analysis types.

## Installation

### Linux (tested on Ubuntu 14.04 LTS and above)

Using `pip`,

`pip install genonets`

In case you get a 'permission' related error, try the following:

`sudo pip install genonets`

You can also install Genonets directly from the source package.

`python setup.py install`

Again, in case you run into permission related errors,

`sudo python setup.py install`

When trying to install genonets on a machine with `Ubuntu 14.04 LTS` that does not already have the required version of 
`python-igraph` installed, `pip` sometimes fails to install the C core of igraph. If that happens, follow these steps:

1. `sudo apt-get install build-essential`
2. `sudo apt-get python-dev`
3. `sudo apt-get install libxml2-dev`
4. `sudo apt-get install libz-dev`
5. `sudo pip uninstall genonets`
6. Finally, `sudo pip install genonets`

### Mac OS X El Capitan

We highly recommend using `virtualenv`, or better yet, `Anaconda`, for installation on Mac OS X El Capitan.

In case you do not already have `virtualenv` installed on your system, use the following command to install 
`virtualenv`:

`pip install virtualenv`

In the directory of your choice, create a virtual environment. In the following example, we will create a virtual 
environment called `venv_genonets`:

`virtualenv venv_genonnets`

Now, activate `venv_genonets` as follows:

`source venv_genonets/bin/activate`

You are now ready to install Genonets. Use the following command:

`pip install genonets`

Note: Every time you need to use `genonets`, you will have to activate the corresponding virtual environment.

### Windows

Instructions for Windows are basically the same, except in certain cases installation of dependencies fails. If that 
happens, follow these steps:

1. Download the 'whl' files for `numpy` and `python-igraph` from http://www.lfd.uci.edu/~gohlke/pythonlibs/. E.g.,
    * `numpy-1.10.2+mkl-cp27-none-win32.whl`
    * `python_igraph-0.7.1.post6-cp27-none-win32.whl`
3. `pip install python_igraph-0.7.1.post6-cp27-none-win32.whl`
4. `pip install numpy-1.10.2+mkl-cp27-none-win32.whl`
5. And finally, `pip install genonets`

## Genonets Server end-of-life

The Genonets Server had to be shut down due to lack of funding for web hosting 
services. Nevertheless, all analyses are still available in the Genonets 
package via this repository. Unfortunately though, the visualization 
features are no longer available. Within the next few weeks, this repository 
will be updated to support Python 3 with significantly improved analyses. Along 
with these changes, the latest documentation will be provided.
