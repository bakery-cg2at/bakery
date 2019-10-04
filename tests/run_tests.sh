#! /bin/sh
#
# run_tests.sh
# Copyright (C) 2016,2017 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#


for d in *; do
    if [ -d $d ] && [ "$d" != "testsuit" ]; then
        cd $d
        ./run_test.sh &> ${d}_log
        RET="$?"
        if [ "$RET" != "0" ]; then
            echo "Error in ${d}"
            cat ${d}_log
            exit $RET
        fi
        cd ..
    fi
done

# Run Python TestCases
cd testsuit
python2 test_preapre_files.py
RET="$?"
cd ..
exit $RET
