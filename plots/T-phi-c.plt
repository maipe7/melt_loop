reset
set terminal qt 0 font 'Calibri,12'
#show terminal

dfile='../outputs/d.txt'

set multi 
set title titlem
set size 0.4,0.9

set xtics nomirror
set ytics nomirror
set grid

unset key
K2C=273.
km2GPa=2800.*10*1e-6
Pa2km=1./2800./10./1e3
maxdepth=80e3

set yrange [80:0]
set ylabel "depth (km)"
set y2label "pressure (GPa)"
set y2tics 0.5
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

imax=30 #30

#test
#pause -1

# load particle files
# particles: <x> <y> <id> <p> <T> <porosity> <peridotite> <peridotiteF> 
file_name(n) = sprintf("../outputs/particles/particles-%05d.0000.gnuplots", n)
array aparticle[imax+1]
do for [iMyr = 0:imax] {aparticle[iMyr+1]=14} # set some value
#aparticle[1]=20; aparticle[11]=20; aparticle[21]=19; aparticle[31]=15 
aparticle[1]=24; aparticle[11]=24; aparticle[21]=23; aparticle[31]=18 
aparticle[1]=25; aparticle[11]=25; aparticle[21]=24; aparticle[31]=19 
#aparticle[1]=22; aparticle[11]=22; aparticle[21]=19; aparticle[31]=16; aparticle[41]=13; aparticle[51]=13
#aparticle[1]=23; aparticle[11]=20; aparticle[21]=13; aparticle[31]=12; aparticle[41]=12 #; aparticle[51]=11
# d.txt: <x> <depth> <T>

## TEMPERATURE ##################################
set origin 0,0
set xrange [0:1200]
set xlabel "temperature (C)"

do for [iMyr = 0:imax:10] {
 p dfile u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l lt 4 lw 1
 #iparticle=5
 #plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every 1:1::iparticle::iparticle w p @PTT ps 1
}
iMyr = 0
 p dfile u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l lt 0 lw 2
 iparticle=aparticle[iMyr+1]
 plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTT ps 1
iMyr = 10
 p dfile u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTT lw 1
 iparticle=aparticle[iMyr+1]
 plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTT ps 2
iMyr = 20
 p dfile u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTT lw 2
 iparticle=aparticle[iMyr+1]
 plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTT ps 3
iMyr = 30
 p dfile u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTT lw 2
 iparticle=aparticle[iMyr+1]
 plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTT ps 4
do for [iMyr = 30:imax:10] {
 p dfile u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTT lw 3
 iparticle=aparticle[iMyr+1]
 plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTT ps 4
}
iMyr=imax
do for [iparticle = 0:100] {
 plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) w p
}

## POROSITY ################################

set origin 0.3,0
set xrange [-0.01:0.5]
set xtics 0.1
set xlabel "melt fraction"

do for [iMyr = 0:imax:10] {
 p dfile u (($4)):(($2)*1e-3) every 1:1::iMyr::iMyr w l lt 4 lw 1
}
iMyr = 0
 p dfile u (($4)):(($2)*1e-3) every 1:1::iMyr::iMyr w l lt 0 lw 2
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTP ps 1
iMyr = 10
 p dfile u ($4):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTP lw 1
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTP ps 2
iMyr = 20
 p dfile u ($4):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTP lw 2
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTP ps 3
iMyr = 30
 p dfile u ($4):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTP lw 2
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTP ps 4
do for [iMyr = 30:imax:10] {
 p dfile u ($4):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTP lw 3
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTP ps 4
}
iMyr=imax
do for [iparticle = 0:100] {
 p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) w lp
}

## COMPOSITION ################################

set origin 0.6,0
set xrange [-0.01:1.5]
set xtics 0.5
set xlabel "composition"

# this is not accurate...
do for [iMyr = 0:imax:10] {
 p dfile u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l lt 4 lw 1
}
iMyr = 0
 p dfile u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l lt 0 lw 2
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTC ps 1
iMyr = 10
 p dfile u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l @LTC lw 1
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTC ps 2
iMyr = 20
 p dfile u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l @LTC lw 2
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTC ps 3
iMyr = 30
 p dfile u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l @LTC lw 2
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTC ps 4
do for [iMyr = 30:imax:10] {
 p dfile u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l @LTC lw 3
 iparticle=aparticle[iMyr+1]
 p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTC ps 4
}
iMyr=imax
do for [iparticle = 0:100] {
 p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) w p
}

unset multi
