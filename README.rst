Genonets
========

This is the Python package used by the `Genonets Server
<http://ieu-genonets.uzh.ch/>`_ for creating and analyzing genotype networks from raw data.

----

Installation
------------

Linux and Mac OS
~~~~~~~~~~~~~~~~

Using pip,

``pip install genonets``

In case you get a 'permission' related error, try the following:

``sudo pip install genonets``

You can also install Genonets directly from the source package.

``python setup.py install``

Again, in case you run into permission related errors,

``sudo python setup.py install``

Windows
~~~~~~~

Instructions for Windows are basically the same, except in certain cases installation of dependencies fails. If that happens, follow these steps:

1. Download the 'whl' files for numpy and igraph from http://www.lfd.uci.edu/~gohlke/pythonlibs/. E.g.,
    i. numpy-1.10.2+mkl-cp27-none-win32.whl 
    ii. python_igraph-0.7.1.post6-cp27-none-win32.whl

3. ``pip install python_igraph-0.7.1.post6-cp27-none-win32.whl``
4. ``pip install numpy-1.10.2+mkl-cp27-none-win32.whl``
5. And finally, ``pip install genonets``

Using Genonets as a command line tool
-------------------------------------

The best way to get started is to work through 'genonets_exmpl_simple.py' available in the 'genonets/sample' directory. The following command can be used to view the list of
command line arguments:

'python genonets_exmpl_simple.py -h'

This directory also includes other sample files, each highlighting
different features.

The details of the analyses used and the attributes computed can be found on the `Learn Genonets
<http://ieu-genonets.uzh.ch/learn>`_ page.