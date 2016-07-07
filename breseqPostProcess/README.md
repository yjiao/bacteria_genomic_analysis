# Post Processing Pipeline for Breseq

Many software packages, such as [Breseq](https://github.com/barricklab/breseq), exist for calling mutations in individual prokaryotic strains. However, jointly calling mutations across sets of related strains is important for detecting artifacts or calling mutations in difficult genomic regions, and is an important part of quality control.

This directory includes scripts used to post-process breseq output from Jiao, Y. and Baym, M. et al, "Population diversity reverses the outcome of antibiotic cycling", 2016, manuscript in progress.

## Usage
**These scripts assume that breseq has been successfully completed on a set of strains.**

```
folder=path/to/folder/containing/all/breseq/output
gbkPath=path/to/reference/bgk/file
./postProcess.sh $folder $gbkPath
```

`postProcess.sh` calls a number of functions:

* `summarize.sh` parses the read count and percent alignment data for all strains.
* `splitEvidence.py` splits read evidence and mutation calls from breseq genome diff (.gd) for easier i/o
* `generateMatrices.py` generates matrices of breseq parameters (such as consensus and polymorphism scores) for called mutations (those that passed all filters). Each row is a mutation, and each column is a strain.
* `forcecall.py` grabs breseq parameters for all called mutations in all strains in the dataset, regardless of whether the mutation passed all filters for any particular strain. Similar to `generateMatrices.py`, `forcecall.py` generates a number of matrices, where each row is a mutation, and each column is a strain.
* `splitCoverageForMatlab.py` and `splitParamsForMatlab.py` reformats data for Matlab.
* `forcecall_allmutations.py` parses information from all "UNKNOWN" locations in breseq. Typically these are areas of very low coverage.



