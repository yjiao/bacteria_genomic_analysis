#!/bin/bash
# This script calculates depth, alignment summaries, and copies relevant files to postProcess outdir

# some jobs take a really long time/ many submissions to complete. This script merges only jobs that are done
# GD tools are buggy and should be avoided when possible. Use intermediate files from direct breseq outputs
#------------------------------------
# Grab file paths
outdir=$1
gbkPath=$2
minReads=10000
minAlign=90
#------------------------------------
# summarize read counts and alignment--this automatically only grabs jobs that have finished
./summarize.sh $outdir
#------------------------------------
mkdir -p $outdir/postProcess/00gd
mkdir -p $outdir/postProcess/01annotated
# concattenate gd file paths that are done
for file in $outdir/*/output/evidence/evidence.gd; do
    oldpath=$file
    strain=${file#*/} #remove first outdir
    strain=${strain%%/*} #remove 2nd to last outdirs
    
    echo "now processing $strain ..."

    newpath=$outdir/$strain/output/$strain.gd
    cp $oldpath $newpath
    cp $newpath $outdir/postProcess/00gd/$strain.gd

    # copy annotation strains
    annotPath=$outdir/$strain/output/evidence/annotated.gd
    annotPath_new=$outdir/postProcess/01annotated/${strain}.gd
    cp $annotPath $annotPath_new
    
    # start more computationally intensive postprocessing scripts
    python splitEvidence.py $outdir # split read evidence and mutation calls for easier i/o
    python mergeAnnot.py $outdir # note this may require manual edits for some files due to strange breseq bug in some cases: copy evidence from 00 to 01 and 02
    
    # Process data for mutations called by breseq
    # This pulls out data for each strain and creates matrices where rows are mutations and columns are strains
    python generateMatrices.py $outdir # generate matrices for all breseq parameters
    python splitCoverageForMatlab.py $outdir/postProcess/04merged/ # reformat coverage data for Matlab
    python splitParamsForMatlab.py $outdir/postProcess/04merged/ # reformat other params for Matlab
    
    # "Force call" data for the set of mutations called in any strain
    python forcecall.py $outdir # grab breseq parameters for all called mutations in all strains in the dataset, output in matrix form
    python splitCoverageForMatlab.py $outdir/postProcess/06forceCalled/ # reformat coverage data for Matlab
    python splitParamsForMatlab.py $outdir/postProcess/06forceCalled/ # reformat other params for Matlab
    python forcecall_allmutations.py $outdir # this grabs all the unknown locations
done
