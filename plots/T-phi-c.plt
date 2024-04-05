reset
set terminal qt 0 font 'Calibri,12'
#show terminal

imax=20
istep=10

ParticlePlot=0 # true/false
ParticlesLast=0 # true/false

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

# load particle files
# particles: <x> <y> <id> <p> <T> <porosity> <peridotite> <peridotiteF> 
file_name(n) = sprintf("../outputs/particles/particles-%05d.0000.gnuplots", n)
imaxx=50
array aparticle[imaxx+1]
do for [iMyr = 0:imaxx] {aparticle[iMyr+1]=14} # set some value
#aparticle[1]=25; aparticle[11]=25; aparticle[21]=24; aparticle[31]=19 
#aparticle[1]=23; aparticle[11]=20; aparticle[21]=13; aparticle[31]=12; aparticle[41]=12 #; aparticle[51]=11
# d.txt: <x> <depth> <T>

## TEMPERATURE ##################################
set origin 0,0
set xrange [0:1200]
set xlabel "temperature (C)"

iMyr = 0
 p dfile u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l lt 0 lw 2
 iparticle=aparticle[iMyr+1]
 if (ParticlePlot) plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTT ps 1
do for [iMyr = istep:imax:istep] {
 isize=iMyr/10
 p dfile u (($3)-K2C):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTT lw isize
 iparticle=aparticle[iMyr+1]
 if (ParticlePlot) {plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTT ps isize}
}
iMyr=imax
if (ParticlesLast) plot file_name(iMyr) u (($5)-K2C):((maxdepth-$2)*1e-3) w p

## POROSITY ################################

set origin 0.3,0
set xrange [-0.01:0.5]
set xtics 0.1
set xlabel "melt fraction"

iMyr = 0
 p dfile u (($4)):(($2)*1e-3) every 1:1::iMyr::iMyr w l lt 0 lw 2
 iparticle=aparticle[iMyr+1]
 if (ParticlePlot) p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTP ps 1
do for [iMyr = istep:imax:istep] {
 isize=iMyr/10
 p dfile u ($4):(($2)*1e-3) every 1:1::iMyr::iMyr w l @LTP lw isize
 iparticle=aparticle[iMyr+1]
 if (ParticlePlot) {p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTP ps isize}
}
iMyr=imax
if (ParticlesLast) p file_name(iMyr) u ($6):((maxdepth-$2)*1e-3) w lp

## COMPOSITION ################################

set origin 0.6,0
set xrange [-0.01:1.5]
set xtics 0.5
set xlabel "composition"

iMyr = 0
 p dfile u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l lt 0 lw 2
 iparticle=aparticle[iMyr+1]
 if (ParticlePlot) p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTC ps 1
do for [iMyr = istep:imax:istep] {
 isize=iMyr/10
 p dfile u (($4*$6+(1.-$4)*$5)):(($2)*1e-3)  every 1:1::iMyr::iMyr w l @LTC lw isize
 iparticle=aparticle[iMyr+1]
 if (ParticlePlot) {p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) every ::iparticle::iparticle w p @PTC ps isize}
}
iMyr=imax
if (ParticlesLast) p file_name(iMyr) u (($6*$8+(1.-$6)*$7)):((maxdepth-$2)*1e-3) w p

unset multi
