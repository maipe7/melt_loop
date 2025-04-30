#reset
set terminal png font 'Calibri,12' size 1000,500
set out model.'.png'

itime=1
shift=0.1

dfile=model.'/ascii_data.txt'

set multi 
set size 0.15,0.9

set xtics nomirror
set ytics nomirror
set grid

unset key
K2C=273.
km2GPa=2800.*10*1e-6
Pa2km=1./2800./10./1e3
maxdepth=80e3

set yrange [0:10]
set ylabel "x (km)"
#set y2label "pressure (GPa)"
#set y2tics 1
#set link y2 via y*km2GPa inverse y/km2GPa
unset ylabel
unset y2label
#set format y ""
#set format y2 ""
set view map
set palette rgbformulae 30,31,32
set cbrange [0:1.1e7]
set colorbox user origin 0,0  size .04,.4
unset colorbox 

set style line 1 lt rgb "black" lw 1 pt 6

# ascii.txt: time           minx           maxx           miny           maxy         volume    temperaturemax_temperature   logviscosity  logviscosityF       porosity   max_porosity             cb         min_cb         max_cb             pc         min_pc         max_pc       velocity   max_velocity      velocityF  max_velocityF    sepvelocitymax_sepvelocity
## TEMPERATURE ##################################
set origin 0,0
set xrange [600:1000]
set xtics 200
set xlabel "temperature"
set title "xxxxxxxxxxxxxxxxxxxxxxxxx   " . model

sp dfile u ($7-K2C):($4/1e3):itime w l palette lw 2

set title " "
unset ytics
## POROSITY ################################

set origin (shift*1),0
set xrange [-0.01:0.3]
set xtics 0.1
set xlabel "melt fraction"

sp dfile u 11:($4/1e3):itime w l palette lw 2

## COMPOSITION ################################

set origin (shift*2),0
set xrange [-0.01:1.5]
set xtics 0.5
set xlabel "composition"

sp dfile u 13:($4/1e3):itime w l palette lw 2

## COMPOSITION SOLID ################################

set origin (shift*3),0
set xrange [-0.01:1.5]
set xtics 0.5
set xlabel "c solid"

sp dfile u 14:($4/1e3):itime w l palette lw 2

## COMPOSITION MELT ################################

set origin (shift*4),0
set xrange [-0.01:10]
set xtics 2
set xlabel "c melt"

sp dfile u 15:($4/1e3):itime w l palette lw 2

## MELT VISCOSITY #############################

set origin (shift*5),0
set xrange [3:7]
set xtics 1
set xlabel "melt viscosity"

sp dfile u 10:($4/1e3):itime w l palette lw 2

## SOLID SHEAR VISCOSITY ##############

set origin (shift*6),0
set xrange [17:21.1]
set xtics 1
set xlabel "shear viscosity"

sp dfile u 9:($4/1e3):itime w l palette lw 2

## DARCY #############################

set origin (shift*7),0
set xrange [-11:-6]
set xtics 1
set xlabel "var Darcy"

sp dfile u (log10($11)*3.0-$10):($4/1e3):itime w l palette lw 2

unset multi

set out 'blank.png'