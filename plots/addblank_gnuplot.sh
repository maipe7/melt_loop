folder=thermal-melt-petrol-6
 awk -f addblanks.awk < outputs/$folder/depth_average.txt > d.txt
 rm -r ./particles/
 mkdir particles
 cp outputs/$folder/particles/*0.0000.gnuplot ./particles/
 gnuplot -p 'T-phi-c.plt'
 #gnuplot -p 'geotherms-T-t.plt'
