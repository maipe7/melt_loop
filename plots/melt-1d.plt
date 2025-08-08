#reset
set terminal png font 'Calibri,12' size 1000,500
set out 'melt-1d.png'

#set term svg; set out "out.svg"; set out 'melt-1d-PTvar.svg'

itime=1
ishift=0
shift=0.1
set size 0.15,0.9
#shift=0.15; set size 0.2,0.9

dfile=model.'/ascii_data.txt'

set multi 

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
set xrange [640:960]
set xtics 200
set xlabel "temperature"
set title "xxxxxxxxxxxxxxxxxxxxxxxxx   " . model

sp dfile u ($7-K2C):($4/1e3):itime w l palette lw 2

set title " "
unset ytics
## POROSITY ################################

ishift = ishift + 1
set origin (shift*ishift),0
set xrange [-0.01:0.3]
set xtics 0.1
set xlabel "melt fraction"

sp dfile u 11:($4/1e3):itime w l palette lw 2

## COMPOSITION ################################

ishift = ishift + 1
set origin (shift*ishift),0
set xrange [-0.01:3.01]
set xtics 0.5
set xlabel "composition"

sp dfile u 13:($4/1e3):itime w l palette lw 2

## COMPOSITION SOLID ################################

ishift = ishift + 1
set origin (shift*ishift),0
set xrange [-0.01:2.1]
set xtics 0.5
set xlabel "c solid"

sp dfile u 14:($4/1e3):itime w l palette lw 2

## COMPOSITION MELT ################################

ishift = ishift + 1
set origin (shift*ishift),0
set xrange [5:15]
set xtics 2
set xlabel "c melt"

sp dfile u 15:($4/1e3):itime w l palette lw 2

## MELT VISCOSITY #############################

ishift = ishift + 1
set origin (shift*ishift),0
set xrange [2.9:5.1]
set xtics 1
set xlabel "melt viscosity"

sp dfile u 10:($4/1e3):itime w l palette lw 2

## log PERMEABILITY #############################
if (1) {
ishift = ishift + 1
set origin (shift*ishift),0
set xrange [-7:-1]
set xtics 1
set xlabel "log10 var permeability"

#sp dfile u (log10($11)*3.0):($4/1e3):itime w l palette lw 2
sp dfile u (log10(($11)**3.0*(1.0-($11))**2)):($4/1e3):itime w l palette lw 2
}

## SOLID SHEAR VISCOSITY ##############
if (0) {
ishift = ishift + 1
set origin (shift*ishift),0
set xrange [16:21.1]
set xtics 1
set xlabel "shear viscosity"

sp dfile u 9:($4/1e3):itime w l palette lw 2
}
## DARCY #############################
if (0) {
ishift = ishift + 1
set origin (shift*ishift),0
set xrange [-11:-6]
set xtics 1
set xlabel "var Darcy"

sp dfile u (log10($11)*3.0-$10):($4/1e3):itime w l palette lw 2
}
unset multi

set out 'blank.png'