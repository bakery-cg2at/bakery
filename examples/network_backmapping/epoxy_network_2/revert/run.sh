#! /bin/sh
#
# run.sh
# Copyright (C) 2015 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#


python start_sim.py --conf hyb_conf.gro --top hyb_topol.top --res_rate 0.0000001 --temperature 300.0 --thermostat_gamma 5.0 --int_step 10 --eq 100 &> sim.log
