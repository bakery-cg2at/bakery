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
        ./run_test.sh
        RET="$?"
        cd ..
        if [ "$RET" != "0" ]; then
            exit $RET
        fi
    fi
done
