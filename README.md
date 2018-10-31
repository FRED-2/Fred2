[![Build Status](https://travis-ci.org/FRED-2/Fred2.svg)](https://travis-ci.org/FRED-2/Fred2) [![Anaconda-Server Badge](https://anaconda.org/bioconda/fred2/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

FRED2 - An Immunoinformatics Framework for Python
-------------------------------------------------
Copyright 2014 by Benjamin Schuber,  Mathias Walzer, Philipp Brachvogel, Andras Szolek, Christopher Mohr, and Oliver Kohlbacher


FRED is a framework for T-cell epitope detection, and vaccine design. It  offers consistent, easy, and simultaneous access to well established prediction methods of computational immunology. FRED can handle polymorphic proteins and offers analysis tools to select, assemble, and design linker sequences for string-of-beads epitope-based vaccines. It is implemented in Python in a modular way and can easily be extended by user defined methods.


Copyright
----------
Fred2 is released under the three clause BSD licence.

Installation
------------

use the following commands:

    $ pip install git+https://github.com/FRED-2/Fred2
    
Dependencies
------------

**Python Packages**
- pandas 
- pyomo>=4.0 
- svmlight 
- MySQL-python>=1.2.4 
- Biopython 
- pyVCF
    
**Thrid-Party Software (not installed through pip)**
   - NetMHC predictor family (NetMHC(pan)-(I/II), NetChop, NetCTL) (http://www.cbs.dtu.dk/services/software.php)
   - PickPocket (http://www.cbs.dtu.dk/services/software.php)
   - Integer Linear Programming Solver (recommended CBC: https://projects.coin-or.org/Cbc)
   
Please pay attanention to the different licensings of the third party tools.

Getting Started
---------------

Users and developers should start by reading our [wiki](https://github.com/FRED-2/Fred2/wiki) and [IPython tutorials](https://github.com/FRED-2/Fred2/tree/master/Fred2/tutorials).
A reference documentation is also available [online](http://fred2.readthedocs.org/en/latest/).

How to Cite
-----------
Please cite:

[Schubert, B., Walzer, M., Brachvogel, H-P., Sozolek, A., Mohr, C., and Kohlbacher, O. (2016). FRED 2 - An Immunoinformatics Framework for Python. Bioinformatics 2016; doi: 10.1093/bioinformatics/btw113](http://bioinformatics.oxfordjournals.org/content/early/2016/02/26/bioinformatics.btw113.short?rss=1)

and the original publications of the used methods.
