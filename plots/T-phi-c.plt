reset
set terminal qt 0 font 'Calibri,12'
#show terminal

set multi 
set size 0.4,0.9

set xtics nomirror
set ytics nomirror
set grid

unset key
K2C=273.
km2GPa=2800.*10*1e-6
Pa2km=1./2800./10./1e3
maxdepth=70e3

set yrange [70:0]
set ylabel "depth (km)"
set y2label "pressure (GPa)"
set y2tics 0.2
set link y2 via y*km2GPa inverse y/km2GPa
unset ylabel
unset y2label
set format y ""
set format y2 ""

PTT='pt 6 lt 8 lw 2'
LTT='lt 8'

PTP='pt 6 lt 7 lw 2'
LTP='lt 7'

PTC='pt 6 lt 6 lw 2'
LTC='lt 6'

#test
#pause -1

# load particle files
imax=30
# particles: <x> <y> <id> <p> <T> <porosity> <peridotite> <peridotiteF> 
file_name(n) = sprintf("particles/particles-%05d.0000.gnuplot", n)
array aparticle[imax+1]
#aparticle[1]=4; aparticle[11]=4; aparticle[21]=3; aparticle[31]=1
#aparticle[1]=5; aparticle[11]=5; aparticle[21]=4; aparticle[31]=2
#aparticle[1]=5; aparticle[11]=5; aparticle[21]=4; aparticle[31]=3
aparticle[1]=10; aparticle[11]=10; aparticle[21]=9; aparticle[31]=5
#aparticle[1]=11; aparticle[11]=11; aparticle[21]=10; aparticle[31]=6
#aparticle[1]=13; aparticle[11]=13; aparticle[21]=13; aparticle[31]=10
# d.txt: <x> <depth> <T>

## TEMPERATURE ##################################
set origin 0,0
set xrange [0:1000]
set xlabel "temperature (C)"
iMyr = 0
 p 'd.txt' u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l lt 0 lw 2
 iparticle=aparticle[iMyr+1]
 plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTT ps 1
iMyr = 10
 p 'd.txt' u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTT lw 1
 iparticle=aparticle[iMyr+1]
 plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTT ps 2
iMyr = 20
 p 'd.txt' u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTT lw 2
 iparticle=aparticle[iMyr+1]
 plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTT ps 3
iMyr = 30
 p 'd.txt' u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTT lw 3
 iparticle=aparticle[iMyr+1]
 plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTT ps 4

## POROSITY ################################

set origin 0.3,0
set xrange [-0.01:0.3]
set xtics 0.1
set xlabel "melt fraction"

iMyr = 0
 p 'd.txt' u (($4)):(($2)*1e-3) every 1:1::iMyr::iMyr w l lt 0 lw 2
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTP ps 1
iMyr = 10
 p 'd.txt' u ($4):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTP lw 1
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTP ps 2
iMyr = 20
 p 'd.txt' u ($4):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTP lw 2
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTP ps 3
iMyr = 30
 p 'd.txt' u ($4):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTP lw 3
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTP ps 4


## COMPOSITION ################################

set origin 0.6,0
set xrange [-0.01:1.4]
set xtics 0.5
set xlabel "composition"

# this is not accurate...
iMyr = 0
 p 'd.txt' u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l lt 0 lw 2
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTC ps 1
iMyr = 10
 p 'd.txt' u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l @LTC lw 1
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTC ps 2
iMyr = 20
 p 'd.txt' u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l @LTC lw 2
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTC ps 3
iMyr = 30
 p 'd.txt' u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l @LTC lw 3
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTC ps 4


unset multi
