# CALCULATION OF SOLIDUS AND LIQUIDUS CURVES,
# MELT FRACTION BASED ON LEVER RULE.
# SINGLE-VALUED COMPOSITION ~ WATER CONTENT (WT%H2O)
# PARAMETRIZATION BASED ON EXPERIMENTS AND THERMOCALC
# FOR FELSIC ROCKS
reset

#set term svg
#set out "out.svg"

# Plotting mode:
# 1=W-T for constant P, one plot
# 2=W-P-T + data
# 3=PT-phi contour plot
# 4=W-T for constant P, W-P for constant T
plotit=4

set tics font "Arial,12"
set xlabel font "Arial, 12"
set ylabel font "Arial, 12"

# liquidus temperature-composition defined as a fit of data for haplogranite melt:
c=600.0
a               = 19. 
b               = 29. 
d               = 360.
e               = 205.
f               = 0.13

# for T=c, w=a*P+b
# for w=0, T-c=d+e*P
# in between, w depends like -T**f 
wliquid(P,T)= (a*P+b)*(1-((T-c)/(d+e*P))**f) #P..GPa,T..C

#fit wliquid(x,y) "w-P-T.dat" using 1:2:3 via a,b,c,d,e,f # xi=0.8 with c=643.
#fit wliquid(x,y) "w-P-T.dat" using 1:2:3 via a,b,d,e,f # xi~=1 with fixed to c=600.

# pressure-dependent wet solidus temperature (Holtz 2001) combined with Thermocalc (Stipska 2019???); works above 0.1 GPa:
Twetsolid(P)=(P<1 ? 640.0 + 150.*(P-1.)**4/1. : 640.0 + 50.*(P-1.)**2/1.) 

# composition where solidus and liquidus curves meet:
#wmax(P)=wliquid(P,Twetsolid(P)) # redundant
#for w=0, dry solidus=dry liquidus - curves cross there:
Tdrysolid(P)= c + d + e*P 

# muscovite dehydration line according to Thermocalc; experiments put it higher (Patino Douce and Harris, 1998):
Tmu(P)=620.0 + 130.*P 

# piecewise linear parametrization of solidus temperature-composition:
# wt%H2O bond in minerals - muscovite, biotite (incl. amphibole..)
#cf. Clemens - pelites 1-1.4 wt%H2O with approx 2x more Bt than Mu, 
#wbt=0.7; wmuJump=0.7; wmu=0.7
# quartzofeldspathic 0.6 wt%H2O with Bt only
#wbt=0.6; wmuJump=0.0; wmu=0.0 # but this is perhaps not realistic - these rocks can crystallize Mu??
# Hasalova, Fig. 12: wbt~0.3, wbt+wmuJump~0.5, wtot~1.0 wt%H2O
#wbt=0.43; wmuJump=0.22; wmu=0.35
# Hasalova, Fig. 9: wbt~0.5, wmuJump~0.2, wmu~0.15
#wbt=0.5; wmuJump=0.2; wmu=0.15
# But Stipska 2019: both Mu and Bt are present - due to high P
#wbt=0.5; wmuJump=0.2; wmu=0.3
# Stipska 2008 granulite?

# quartzofeldspathic:
wbt=0.5; wmuJump=0.2; wmu=0.3
Wtot = 0.9 # wtot: stipska 0.7, 1.0; hasalova 0.7-0.8; guilmette 0.5, 1.2

# pelitic:?
#wbt=0.75; wmuJump=0.3; wmu=0.45
#Wtot = 1.35

# total water content in saturated rock:
wsat=wmu+wmuJump+wbt
pr 'wsat=',wsat

# temperature interval over which wet solidus and muscovite dehydration happen, respectively:
DT=10.

# reference pressure (used only when P-dependence of wbt,wmu,wmuJump is assumed):
P0=1. # GPa

wlin0(P,T)=(wbt)*(T-Tdrysolid(P))/(Tmu(P)-Tdrysolid(P)) # fixed (zeroth approximation)
#wlin0(P,T)=(wbt)*(T-Tdrysolid(P))/(Tmu(P0)-Tdrysolid(P)) # same reference point - preferred
#wlin0(P,T)=(wbt)*(T-Tdrysolid(P))/(Tmu(P0)-Tdrysolid(P0)) # fixed slope
wlin1(P,T)=(wlin0(P,Tmu(P))+wmuJump)+(wmu)*(T-Tmu(P))/(Twetsolid(P)-Tmu(P)) # fixed starting point (zeroth approximation)
#wlin1(P,T)=(wlin0(P,Tmu(P))+wmuJump)+(wsat-wmuJump-wlin0(P,Tmu(P)))*(T-Tmu(P))/(Twetsolid(P)-Tmu(P)) # fixed saturation point
#wlin1(P,T)=(wlin0(P,Tmu(P))+wmuJump)+(wsat-wmuJump-wbt)*(T-Tmu(P))/(Twetsolid(P0)-Tmu(P0)) # fixed slope - preferred

wsolid(P,T)=( T > Twetsolid(P)-DT ? \
  ( T > Twetsolid(P) ? \
  ( T > Tmu(P)-DT ? \
  (T>Tmu(P) ? \
   wlin0(P,T) : \
   wlin0(P,Tmu(P))+wmuJump*(Tmu(P)-T)/DT ) : \
   wlin1(P,T) ): \
  wlin1(P,Twetsolid(P))+(wliquid(P,Twetsolid(P)-DT)-wlin1(P,Twetsolid(P)))*((Twetsolid(P)-T)/DT)) : \
  wliquid(P,T))

#################################
###### PLOTTING #################
#################################

if (plotit==1) {
pr plotit
set title "Water content at solidus and liquidus temperature"
set samples 180
set isosamples 180

set parametric
set trange [600:1200]
#set trange [Tmin(1.0):1300]

set ylabel "T(C)"
set yrange [600:1200]
set xlabel "wt%H2O"
set xrange [0:20]
p wliquid(1.3,t),t lw 3 lt 1 dt 2 t "P=1.5 GPa, liquidus", wsolid(1.5,t),t lw 3 lt 1 t "solidus", \
 wliquid(1.0,t),t lw 3 lt 2 dt 2 t "P=1.0 GPa, liquidus", wsolid(1.0,t),t lw 3 lt 2 t "solidus", \
 wliquid(0.5,t),t lw 3 lt 3 dt 2 t "P=0.5 GPa, liquidus", wsolid(0.5,t),t lw 3 lt 3 t "solidus"

unset parametric
}
if (plotit==2) {
#################################
# solidus and liquidus in T-P-W space + data for liquidus:

set xrange [1.5:0.0]
set yrange [600:1300]
set zrange [0:25]
set xlabel "pressure (GPa)"
set ylabel "temperature (C)"
set zlabel "wt%H2O"
#set contour
set cntrparam levels incremental 1,1,20
set samples 80
set isosamples 80

sp wliquid(x,y), "w-P-T.dat" w p pt 7, wsolid(x,y)
}
#################################
# Melt fraction (P,T) contour plot for fixed water content
if(plotit==3){
phi(P,T,W)=(wsolid(P,T)-W)/(wsolid(P,T)-wliquid(P,T))
set multi

set xlabel "Temperature (C)"
set xrange [600:1000]
set ylabel "Pressure (GPa)"
set yrange [0:2]
set view map
unset surface
set pm3d
set palette defined (0 0.9 0.9 0.9, 1 0.9 0.9 0.9)
set zrange [0:1.0]
#set cbrange [0:1.0]
set contour
unset colorbox
set cntrparam levels discrete 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5
set style textbox opaque fc "#00EEEEEE" noborder

set origin 0,0; set margins screen 0.2, screen 0.7, screen 0.2, screen 0.7; set nokey

Title = sprintf("%s %2.1f %s", "Melt fraction for fixed", Wtot,"wt%H2O")
set title Title
set samples 200; set isosamples 200

# create internal table with contours:
set table $Contour
    splot phi(y,x,Wtot)
unset table

splot phi(y,x,Wtot)
unset pm3d
unset tics; unset xlabel; unset ylabel; unset title

plot $Contour w l lt -1 lw 3
plot $Contour u 1:2:3 every 200::70 w labels boxed notitle 
#plot $Contour u 1:2:3 every 200::70 w labels font "Arial,12" boxed notitle 

unset multi

}

#################################
# solidus and liquidus lines: W-T for constant P, W-P for constan T:
if(plotit==4){
set multi layout 2,3

set parametric
set trange [600:1200]
#set trange [Tmin(1.0):1200]
unset title

set ylabel "T(C)"
set yrange [600:1200]
#set yrange [Tmin(1.0):1300]
set xlabel "wt%H2O"
set xrange [0:10]
P=0.5
p wliquid(P,t),t lw 3 lt 1 dt 2 t "P=0.5 GPa, liquidus", wsolid(P,t),t lw 3 lt 1 t "solidus"
P=1.0
p wliquid(P,t),t lw 3 lt 2 dt 2 t "P=1.0 GPa, liquidus", wsolid(P,t),t lw 3 lt 2 t "solidus"
P=2.0
p wliquid(P,t),t lw 3 lt 3 dt 2 t "P=2.0 GPa, liquidus", wsolid(P,t),t lw 3 lt 3 t "solidus"

set trange [0:2.0]
set ylabel "P(GPa)"
set yrange [0:2.0]
set xlabel "wt%H2O"
set xrange [0:10]
T=650.
p wliquid(t,T),t lw 3 lt 1 dt 2 t "T=650 C, liquidus", wsolid(t,T),t lw 3 lt 1 t "solidus"
T=700.
p wliquid(t,T),t lw 3 lt 2 dt 2 t "T=700 C, liquidus", wsolid(t,T),t lw 3 lt 2 t "solidus"
T=800.
p wliquid(t,T),t lw 3 lt 3 dt 2 t "T=800 C, liquidus", wsolid(t,T),t lw 3 lt 3 t "solidus"
unset parametric
unset multi
}

#set out
set term win

##############################