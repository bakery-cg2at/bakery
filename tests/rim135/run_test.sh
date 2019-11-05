#! /bin/sh
#
# run_test.sh
# Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#


python ../../src/prepare_files.py --options settings.xml --verbose

# Compare files and check if they are okey
diff ref_hyb_topol.top hyb_topol.top
TOPOL=$?
diff ref_hyb_conf.gro hyb_conf.gro
CONF=$?

[ "$TOPOL" = "0" ] && [ "$CONF" = "0" ] && echo "$0 OK" || echo "$0 Fail"

if [ "$TOPOL" = "0" ] && [ "$CONF" = "0" ]; then
    rm -f hyb_conf.gro hyb_topol.top
    rm -f missing_definitions.txt
    exit 0
else
    exit 1
fi
