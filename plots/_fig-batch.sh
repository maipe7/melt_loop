#!/bin/bash

#date=$(date +%Y%m%d)
#outdir="figs"$date
#mkdir $outdir

for file in thermal-melt-petrol-2d-varCPS1-30kmdeep-200-10-k8-eta1-P* #box-eta-ccx10-3-r125* #[0-2]*
#for file in ./oF0* #[0-2]*
do
 pvpython ./postproc-melt-FigM.py $file
done

