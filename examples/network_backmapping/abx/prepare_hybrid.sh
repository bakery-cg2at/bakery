#!/bin/bash

CR_DIR=$1

for t in $CR_DIR/*; do
    if [ -f $t ]; then
        continue
    fi
    oldpwd=$PWD
    cd $t
    dname=`basename $t`
    mkdir -v -p data
    prepare_files.py --options backmapping_settings.xml --generate-only-graph
    python2 ../../scan_bonds.py 
    prepare_files.py --options backmapping_settings.xml

    qsub sub.pbs -N reverse_${dname}

    cd $oldpwd
done
