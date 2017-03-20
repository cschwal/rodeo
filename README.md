# RODEO
Rapid ORF Description and Elucidation Online

For usage of RODEO, see [http://www.ripprodeo.org](http://www.ripprodeo.org)

## Contents
The present directory contains the following files related to RODEO and its uses:

### Core RODEO scripts
- rodeo.pl: the main RODEO Perl script itself
- toPeptides.pl: helper script used by RODEO
- translate.pl: helper script used by RODEO

### Diagnostic
- installation_test.sh: shell script to test if RODEO is properly installed

### Machine learning scripts
- svmoptimize.py: Python script for use with scikit-learn in order to automate parameter optimization
- svmclassify.py: Python script for use with scikit-learn in order to learn from one data set and classify another

### Related utilities
- RODEOmulti.sh: shell script to successively run RODEO on all .inp files within a folder
- giToAcc.pl: Perl script to convert a list of GI numbers into GenBank accession identifiers

#### Notes
- RODEO is maintained by the [Mitchell Research Lab](http://www.scs.illinois.edu/mitchell/)
