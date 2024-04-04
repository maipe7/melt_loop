folder=../outputs/thermal-melt-petrol-chwide-bm-k0-r2
#folder=thermal-melt-petrol-chwide-th-et11-c06-cool-k7
 awk -f addblanks.awk < $folder/depth_average.txt > ../outputs/d.txt
 rm -r ../outputs/particles/; mkdir ../outputs/particles
 #cp $folder/particles/*.0000.gnuplot ./particles/
pwd
 cp $folder/particles/*0.0000.gnuplot ../outputs/particles/

for particlefile in ../outputs/particles/*; do
 sort -t' ' -b -nk3 $particlefile --output=${particlefile}s
done
 gnuplot -e "titlem='${folder: 31}'" -p 'T-phi-c.plt'
 