#!/bin/bash

cd ../outputs
mkdir plots-all

for f in *; do
  #if [ -d "$f" ]; then
  if [[ $f == *"melt-1d"* ]] && [ -d "$f" ]; then
    echo $f
    echo "reset; model='${f}'" > melt-1d-all.plt
    cat ../plots/melt-1d.plt >> melt-1d-all.plt  
    gnuplot melt-1d-all.plt
    cp melt-1d.png plots-all/${f}.png
  fi
done


