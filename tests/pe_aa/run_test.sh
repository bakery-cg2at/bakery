#! /bin/sh
#
# run_test.sh
# Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#


python ../../src/prepare_files.py --options backmapping_pe_aa.xml

# Compare files and check if they are okey
diff -q ref_hyb_topol.top hyb_topol.top
TOPOL=$?
diff -q ref_hyb_conf.gro hyb_conf.gro
CONF=$?

if [ "$TOPOL" = "0" ] && [ "$CONF" = "0" ]; then
    rm -f hyb_conf.gro hyb_topol.top
fi

[ "$TOPOL" = "0" ] && [ "$CONF" = "0" ] && echo "$0 OK" || echo "$0 Fail"
[ "$TOPOL" = "0" ] && [ "$CONF" = "0" ] && exit 0 || exit 1
