#!/bin/bash

CR=$1

for d in $CR/p7abx_*; do
    if [ ! -d $d ]; then
        continue
    fi
    echo $d
    dname=`basename $d`
    eqdir=${d}/data/eq
    mkdir -p $eqdir
    cp -rvf template_eq/* $eqdir
    cp -v ${d}/data/*topol_final_aa*.top $eqdir/topol.top
    cp -v ${d}/data/*final_aa_0.0005*.gro $eqdir/conf.gro
    sed -i 's/.*cg_ffnb.itp.*//g' $eqdir/topol.top
    oldpwd=$PWD
    cd $eqdir
    qsub -N eq_${dname} MyJobScript.sh
    cd $oldpwd
done
