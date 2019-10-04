#!/bin/bash
set -euo pipefail

CR_DIR=$1

for t in $CR_DIR/p7abx*.top; do
    python2 convert_topology.py $t `dirname $t`/cg_`basename $t`
done
    
for t in $CR_DIR/cg_*.top; do
    dname=`dirname $t`
    fname=`basename $t .top`
    bname=`echo $fname | cut -f2-3 -d_`
    echo $dname $fname $bname
    dsdir=$CR_DIR/$bname
    mkdir -v -p $dsdir

    cp -v $t $dsdir/cg_topol.top
    cp -v $dname/${bname}_confout.gro $dsdir/cg_conf.gro
    cp -v ./backmapping_settings.xml $dsdir/
    cp -v ./{m1.gro,m1.itp,cg_ffnb.itp,params} $dsdir/
    cp -v ./*.gro $dsdir/
    cp -v ./*.itp  $dsdir/
    cp -v ./sub.pbs $dsdir/

    for tab in ./table_*; do
        ln -s `pwd`/$tab $dsdir/
    done
done


