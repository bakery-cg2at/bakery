#!/bin/bash

DS=$1
REACTION_DATA=~/staging/abx/reactions/ds1_new/$DS


for f in $REACTION_DATA/p7abx_*[0-9]_output_topol.top; do
#for f in $REACTION_DATA/p7abx_3435_output_topol.top; do
    cp -v $f $DS
done

for f in $REACTION_DATA/p7abx_*[0-9]_confout.gro; do
#for f in $REACTION_DATA/p7abx_3435_confout.gro; do
    cp -v $f $DS
done
