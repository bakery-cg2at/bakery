#!/bin/bash

CR=$1

for d in $CR/p7abx_*; do
    if [ ! -d $d ]; then
        continue
    fi
    echo $d
    dname=`basename $d`
    eqdir=${d}/data/eq/nvt
    if [ ! -d $eqdir ]; then
        continue
    fi
    echo $eqdir $dname
    mkdir -p results/$CR
    cp -v $eqdir/confout.gro results/$CR/${dname}_confout.gro
    cp -v $eqdir/topol.top results/$CR/${dname}_topol.top
done
