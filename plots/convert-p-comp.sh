#!/bin/bash

pfolder=.
folder=thermal-melt-petrol-2d-varCPS1-30kmdeep-200-10-k8-eta1-
#convert scales-p-.png scales-c-.png +append scales-pc.png

 #folder=${pfolder:2}
 echo $folder
 for i in $(seq -f "%04g" 0 80)
 do
   echo $i
   convert porosity/$folder/$folder.$i.png composition/$folder/$folder.$i.png +append pc/$folder.$i.png
   convert pc/$folder.$i.png scales-pc-.png -append pc-scales/$folder.$i.png
   ##convert p/$folder/$folder.$i.png u/$folder/$folder.$i.png +append pu/$folder.$i.png
   ##convert pu/$folder.$i.png c/$folder/$folder.$i.png +append puc/$folder.$i.png
   ##convert puc/$folder.$i.png scales-white.png -append puc-scales/$folder.$i.png
 done

# ffmpeg -r 12 -i $pfolder/pc/$folder.%04d.png $folder.wmv
# ffmpeg -r 12 -i $pfolder/pc/$folder.%04d.png -c:v copy $folder.avi
 ffmpeg -framerate 12 -i $pfolder/pc-scales/$folder.%04d.png  -c:v libx264 -crf 0 $folder.mp4
