folder0=../outputs
folder1=thermal-melt-bm-k8-v0
folder=$folder0/$folder1
 awk -f addblanks.awk < $folder/depth_average.txt > $folder0/d.txt
 rm -r $folder0/particles/; mkdir $folder0/particles
 cp $folder/particles/*.0000.gnuplot $folder0/particles/
 #cp $folder/particles/*0.0000.gnuplot $folder0/particles/

for particlefile in $folder0/particles/*; do
 sort -t' ' -b -nk3 $particlefile --output=${particlefile}s
done
 gnuplot -e "titlem='${folder: 31}'" -p 'T-phi-c.plt'
 