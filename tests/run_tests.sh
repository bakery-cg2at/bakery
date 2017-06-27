#! /bin/sh
#
# run_tests.sh
# Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#


for d in *; do
    if [ -d $d ]; then
        cd $d
        ./run_test.sh | tee log
        RET="$?"
        if [ "$RET" != "0" ]; then
            echo "Error"
            cat log
            exit $RET
        fi
        cd ..
    fi
done
python testsuit/test_preapre_files.py
