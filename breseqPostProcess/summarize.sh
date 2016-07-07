#!/bin/bash
outdir=$1
outdir=${outdir////}
grep 'trimmed-pair1</td><td align="right">' $outdir/*/output/summary.html | while read line; do
    strain=${line%%:*}
    strain=${strain#*/}
    strain=${strain%%/*}

    reads=${line#*right\">}
    reads=${reads%%\<*}
    reads=${reads//,/}
    
    align=${line#*bases}
    align=${align#*bases}
    align=${align#*right\">}
    align=${align%%\%*}
    
    echo -e "${strain}\t${reads}\t${align}"
done > $outdir/postProcess/summary.txt

