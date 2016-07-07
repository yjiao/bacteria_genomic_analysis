# Genomic Analysis for Bacteria

This repository provides scripts used to extend the functionality of [Breseq](https://github.com/barricklab/breseq), a popular package for calling mutations in prokaryotic organisms.

* **breseqPostProcess** contains scripts used to parse and collate breseq output data for large datasets containing many strains (hundreds to thousands). It additionally pulls out all breseq parameters (such as polymorphism and consensus scores) for every variant called in at least one strain over the entire dataset, which is important for downstream quality control.

* **samtoolsForceCall** Similar to **breseqPostProcess**, scripts in this directory takes a list mutations and a list of bam files, and pulls out alternate read count and frequency for every mutation in every bam file. Scripts in this directory can be easily adapted to other pipelines outside of breseq.

* **matlabGUI** Downstream processing and visualization in Matlab. Unlike the previous two directories, scripts in this directory are project-specific and may take extensive work for other projects. However, it also includes a nice **interactive GUI** which takes the output files of **breseqPostProcess**, fully integrating the outputs of breseq, and allows users to visualize relationships between different breseq parameters in their sequencing data.

A screenshot of the breseq viewer. Breseq parameters are displayed on the right part of the window. Any parameter can be scattered with respect to another by choosing it as either the x or y variable. Parameters with the "FC" suffix indicate these parameters are of variants that did not pass all breseq filters, but were parsed from **breseqPostProcess**.

The plot on the right displays "FC" parameters as gray dots (variants not passing all filters), and "real" variants (passing all filters) as blue dots.

![Viewer](viewer1.png)

Clicking on any blue dot brings up the breseq HTML page, which displays information on the mutation as well as read alignments.
![Viewer](viewer2.png)

Clicking on a particular point will also display relevant data on the lower right corner of the window.
![Viewer](viewer4.png)

Choosing the same variable for the x and y coordinates draw histograms for the parameter values.
![Viewer](viewer3.png)

