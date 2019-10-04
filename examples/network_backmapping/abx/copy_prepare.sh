#!/usr/bin/bash

python2 convert_topology.py reaction/p7abx_*_output_topol.top cg_topol.top
cp reaction/p7abx_*[0-9]_confout.gro cg_conf.gro
prepare_files.py --options backmapping_settings.xml --allow-no-bonds
rm _*
