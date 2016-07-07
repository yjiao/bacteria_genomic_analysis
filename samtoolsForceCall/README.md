# Force Call All Substitutions

Many software packages, such as [Breseq](https://github.com/barricklab/breseq), exist for calling mutations in individual prokaryotic strains. However, many experiments, such as multi-step evolutionary experiments, require jointly calling mutations across sets of related strains. This is important for detecting artifacts or calling mutations in difficult genomic regions. Particularly in evolutionary experiments, given that a variant exists in a later population, one may want to find out whether this variant as present at low levels in the ancestor.

These scripts provide a way to measure the number of reads supporting specific mutations in a given list. Given a list of bam files from a set of strains, and a list of mutations, these scripts return the count and frequency of alternative alleles of each mutation in each bam file.

The current scripts were written with the directory structure of breseq in mind. However, it can be easily adapted for other sequencing pipelines.

**NOTE** The current version is a simple tool for counting substitutions only. It cannot call insertions and deletions (indels) at the moment.


## Usage
**mutationlist.csv** should be a comma delimited file of chromosomal position, reference, and alternate alleles.

Sample mutationlist.csv:

```
2987,G,A
5221,A,G
5286,G,A
5287,G,C
5299,A,G
5300,T,G
5306,A,T
5338,T,G
```

To force call all mutations in all strains:

```
breseqFolder=path/to/breseq/output
refpath=path/to/ref/fasta/file
mutationlist=path/to/list/of/mutations
./runForceCall.sh $breseqFolder $refpath $mutationList
```
This outputs two files:

* forcecalled_counts.txt
* forcecalled_frequencies.txt

For each file, a row corresponds to a single **mutation** in the mutationlist file. Each column corresponds to a single **strain** in the strainlist.

## History
2016-07-07 Initial commit
