## PLOTTING ROUTINE

# this version of the original script (plot_simulation_experiment_mult_comparison.plt) receives the date(i) arguments in form 
# of expi from the command line; also the modei is passed

# PLOTS EXPERIMENTS LISTED IN date(i) ACCORDING TO THE MODES IN modes(i)
# PLOTS SIMULATIONS IN THE FOLDERS 'simpath/experiment_0N' where N is the integer 1,2,3, etc
# WHERE THE STRING simpath CAN BE DEFINED SO AS TO POINT TO THE DESIRED LOCATION
# PLOTS OPTIONALLY ALSO THE SIMULATIONS LOCATED IN THE FOLDERS 'experiments_0N<comflag>' WHERE
# THE STRING 'comflag' CAN BE SPECIFIED TO POINT TO THE LOCATION
# PLOTS OPTIONALLY ALSO THE EQUILIBRIUM THEORY SOLUTION W VELOCITY VARIATIONS IF THE FLAG 
# eqflag is set to 1.

simpath = ''
# nexp = 2	


comflag = ''
eqflag = 0 # 1 or 0
rhplot = 0


ymax = 10
tmaxh = 10
date(i) = i == 1 ? exp1 : \
		  i == 2 ? exp2 : \
		  i == 3 ? exp3 : \
		  i == 4 ? exp4 : \
		  i == 5 ? exp5 : \
		  i == 6 ? exp6 : \
		  i == 7 ? exp7 : \
		  exp8
 
modes(i) = i == 1 ? mode1 : \
		   i == 2 ? mode2 : \
		   i == 3 ? mode3 : \
		   i == 4 ? mode4 : \
		   i == 5 ? mode5 : \
		   i == 6 ? mode6 : \
		   i == 7 ? mode7 : \
		   mode8
		   
color(n) = n == 1 ? "blue" : \
		   n == 2 ? "red"  : \
		   n == 3 ? "steelblue" : \
		   n == 4 ? "dark-red" : \
		   n == 5 ? "web-green" : \
		   n == 6 ? "bisque" : \
		   n == 7 ? "aquamarine" : \
		   "dark-olivegreen"

destpath = '\\d.ethz.ch\dfs\groups\mavt\spl\Temp\adsorption\SmallColumn\eval\experiments\experiments_by_date\'
destgnu(i) = sprintf('%s%s\',destpath,date(i))
fileid(i) = sprintf('%s_%s',date(i),modes(i))
# MS 
fMSname(i) = sprintf('%soutputMS_%s.dat',destgnu(i),fileid(i))
# Labview, x = 'means', 'Tcol', 'column'
fLVnamefunc(i,x) = sprintf('%soutputLV_%s_%s.dat',destgnu(i),x,fileid(i))
# pressure drop
fdPname(i) = sprintf('%sdP_%s.dat',destgnu(i),fileid(i))

# simulation details
paths(i,inp) = sprintf('experiment_0%i%s/exitprofile.txt',i,inp)
pathl(i,inp) = sprintf('experiment_%i%s/exitprofile.txt',i,inp)
path(i,inp) = i<10?paths(i,inp) : pathl(i,inp)
pathtl(i,inp) = sprintf('experiment_0%i%s/temperatures.txt',i,inp)
pathth(i,inp) = sprintf('experiment_%i%s/temperatures.txt',i,inp)
patht(i,inp) = i<10?pathtl(i,inp) : pathth(i,inp)

# equilibrium theory with velocity (EQV)
eqvelpath = 'D:\Users\mhefti\Documents\Projects\dynamic_experiments_modelling\equilibrium_theory\eq_theory_velocity_v00_experiments\output'
path_eq_vel(i) = sprintf('%s\Exp_%s_adsorption.dat',eqvelpath,i)

# m x n matrix of plots 
m = 2
n = 1

# Dimensions
Awidth = 18
Aheight = Awidth*sqrt(2)


# rh function
a1 =- 7.85951783; a2 =  1.84408259; a3 =-11.7866497; a4 = 22.6807411; a5 =-15.9618719; 
a6 =  1.80122502; 
Tcrit = 647.096			# K critical temperature
Pcrit = 22.064*10**6	# bar critical pressure
Pvap(T)=Pcrit*exp((Tcrit/T)*(a1*(1-(T/Tcrit))+a2*(1-(T/Tcrit))**1.5+a3*(1-(T/Tcrit))**3+a4*(1-(T/Tcrit))**3.5+a5*(1-(T/Tcrit))**4+a6*(1-(T/Tcrit))**7.5))

set macros

############################################################################################################################################
# MACROS 																							    
############################################################################################################################################

# labeling
YYLABEL = "set ylabel 'mole fraction [%]'"
TYLABEL = "unset ylabel; \
		   set ylabel 'temperature [°C]'"

YXLABEL = ""
TXLABEL = "unset xlabel; \
		   set xlabel 'time [hours]"

# environment		   
tsize = .4
MARGTIC  = "set border lw 4; \
			set tics scale tsize; \
			set grid"		
						
L0 = 0.02
B0 = 0.0

pwidth = (1.- L0)/n
pheight = (1. - B0)/m
pheightt = pheight*.65
xpos(in) = L0 + (in-1)*pwidth;
yposy    = B0 + pheight;
ypost =  yposy - 1.08*pheightt + L0

			
YYTICS   = "set ytics autofreq scale tsize format '%.1f'"
TYTICS   = "set ytics autofreq 1 scale tsize format '%.1f'"

YXTICS   = "set xtics autofreq scale tsize format ''"

Tmaxh = 0
Tminh = 0

do for [i=1:nexp] {
	stats fLVnamefunc(i,'Tcol') u 1:2 nooutput
	# get the maximum time
	# tmaxn = STATS_max_x/3600
	# if (tmaxn > tmaxh) {
	#	tmaxh = tmaxn } 

	# get the maximum temperature
	Tmaxn = STATS_max_y
	Tminn = STATS_min_y
	
	if (Tmaxn > Tmaxh) {
		Tmaxh = Tmaxn } 
	
	if (Tminn > Tminh) {
		Tminh = Tminn }
		
}

Tminh = Tminh - 1
Tmaxh = Tmaxh + 1



stats fLVnamefunc(1,'Column') u 1:8 nooutput
skiplines = STATS_index_min_y
skiplines = 1

if (tmaxh < 10) {
	TXTICS   = "set xtics autofreq scale tsize format '%.1f'"
} else {
	TXTICS   = "set xtics autofreq scale tsize format '%.0f'"
}

# plotting experiments

if (rhplot == 0) { # plot against mole fraction

YYEXIT1 = "every ::skiplines u ($1/3600):($12*100) title 'experimental' axes x1y1 w l lt -1 lw 6"
YYEXITEQ1VEL = "u ($1/3600):($2) title 'EQT-V' axes x1y1 w l lt 17 lw 6"
YCEXIT1 = "every ::skiplines u ($1/3600):($14*100) title '' axes x1y1 w l lt 0 lw 6"
YLVEXIT = "every ::skiplines u ($1/3600):($8*100) title 'experimental' axes x1y1 w l lt 1 lw 6"
TMID1   = "every ::skiplines u ($1/3600):($2) axes x1y1 w l lt -1 lw 6"

# plotting simulations
YYEXIT1SIM = "u (($1 + skiplines)/3600):($2*100) title 'detailed' axes x1y1 w l lt -1 lw 6"
YYEXIT1SIMCOMP = "u (($1 + skiplines)/3600):($2*100) title '' axes x1y1 w l lt 2 lw 6"
} else {  # plot as relative humidity

## merge functions
YYEXIT1 = "every ::skiplines u ($1/3600):($12*100) title 'experimental' axes x1y1 w l lt -1 lw 6"
YYEXITEQ1VEL = "u ($1/3600):($2) title 'EQT-V' axes x1y1 w l lt 17 lw 6"
YCEXIT1 = "every ::skiplines u ($1/3600):($14*100) title '' axes x1y1 w l lt 0 lw 6"
YLVEXIT = "every ::skiplines u ($1/3600):(100.0*1e5*$8*$2/Pvap($6)) title 'experimental' axes x1y1 w l lt 1 lw 6"
TMID1   = "every ::skiplines u ($1/3600):($2) axes x1y1 w l lt -1 lw 6"

# plotting simulations
# combine 
YYEXIT1SIM = "u (($1 + skiplines)/3600):($2*$8*100/Pvap($11)) title 'detailed' axes x1y1 w l lt -1 lw 6"
#YYEXIT1SIM = "u (($1 + skiplines)/3600):($2*100) title 'detailed' axes x1y1 w l lt -1 lw 6"
YYEXIT1SIMCOMP = "u (($1 + skiplines)/3600):($2*100) title '' axes x1y1 w l lt 2 lw 6"
}



TMID1SIM  = "u (($1 + skiplines)/3600):($2-273.15) title '' axes x1y1 w l lt -1 lw 6"
TMID1SIMCOMP  = "u (($1 + skiplines)/3600):($2-273.15) title '' axes x1y1 w l lt 2 lw 6"

# size of plots
YSIZE   = "set size pwidth,pheight"
TSIZE   = "set size pwidth,pheightt"

set terminal pdfcairo noenhanced dashed color font "CMU Sans Serif,16" size Awidth cm, Aheight cm

outfile = simpath.'experiments_comparison_'.comflag.'.pdf'
set output outfile

# extract current date:
rhdate = 140714

set style fill transparent solid 0.5 border
#----------------------------------------------------------------------------------------------------
set multiplot layout n,m rowsfirst
#----------------------------------------------------------------------------------------------------

set title '' offset 0,-.5
set key noauto
set grid

#####################################################################################################
# Figure 1,1
#####################################################################################################
#@GFX1Y
ifig = 1
jfig = 1
set origin xpos(ifig),yposy

@MARGTIC
@YYLABEL
@YXLABEL
@YYTICS
@YXTICS
@YSIZE

set yrange [0:*]
set xrange [0:*]


do for [i=1:nexp] {
	set key inside right bottom
	set xrange [0:tmaxh]
	set yrange [0:ymax]
	currdate = date(i)
	set origin xpos(ifig),yposy
	@YSIZE

	# 'experiments_by_date/140730/adsorption/'
	if (rhplot == 0) {
		if (currdate >= rhdate && eqflag == 0){
			plot fLVnamefunc(i,'column')	@YLVEXIT lc rgb color(i),\
				simpath.path(i,'')			@YYEXIT1SIM lc rgb color(i),\
				simpath.path(i,comflag)	@YYEXIT1SIMCOMP lc rgb color(i)
			}
		if (currdate >= rhdate && eqflag == 1){
			plot fLVnamefunc(i,'column')	@YLVEXIT lc rgb color(i),\
				simpath.path(i,'')			@YYEXIT1SIM lc rgb color(i),\
				simpath.path(i,comflag)	@YYEXIT1SIMCOMP lc rgb color(i),\
				path_eq_vel(currdate)		@YYEXITEQ1VEL lc rgb color(i)
			}
		if (currdate < rhdate && eqflag == 0){
			plot fMSnamefunc(i)				@YYEXIT1 lc rgb color(i),\
				fMSnamefunc(i)				@YCEXIT1 lc rgb color(i),\
				simpath.path(i,'')			@YYEXIT1SIM lc rgb color(i),\
				simpath.path(i,comflag)	@YYEXIT1SIMCOMP lc rgb color(i)
			}
		if (currdate < rhdate && eqflag == 1){
			plot fMSnamefunc(i)				@YYEXIT1 lc rgb color(i),\
				fMSnamefunc(i)				@YCEXIT1 lc rgb color(i),\
				simpath.path(i,'')			@YYEXIT1SIM lc rgb color(i),\
				simpath.path(i,comflag)	@YYEXIT1SIMCOMP lc rgb color(i),\
				path_eq_vel(currdate)		@YYEXITEQ1VEL lc rgb color(i)
			}
	} else{
		if (currdate >= rhdate && eqflag == 0){
			fname1 = simpath.path(i,'')
			fname2 = simpath.patht(i,'')
			mergefile = simpath.'mergefile.dat'
			syscommand = sprintf('%s.pl %s %s %s','concat',fname1,fname2,mergefile)
			system syscommand
			#system 'horzcat2.pl experiment_01/exitprofile.txt experiment_01/temperatures.txt mergefile.dat'
			plot fLVnamefunc(i,'column')	@YLVEXIT lc rgb color(i),\
				 mergefile @YYEXIT1SIM lc rgb color(i),\
				simpath.path(i,comflag)	@YYEXIT1SIMCOMP lc rgb color(i)
			}
		if (currdate >= rhdate && eqflag == 1){
			plot fLVnamefunc(i,'column')	@YLVEXIT lc rgb color(i),\
				simpath.path(i,'')			@YYEXIT1SIM lc rgb color(i),\
				simpath.path(i,comflag)	@YYEXIT1SIMCOMP lc rgb color(i),\
				path_eq_vel(currdate)		@YYEXITEQ1VEL lc rgb color(i)
			}
		if (currdate < rhdate && eqflag == 0){
			plot fMSnamefunc(i)				@YYEXIT1 lc rgb color(i),\
				fMSnamefunc(i)				@YCEXIT1 lc rgb color(i),\
				simpath.path(i,'')			@YYEXIT1SIM lc rgb color(i),\
				simpath.path(i,comflag)	@YYEXIT1SIMCOMP lc rgb color(i)
			}
		if (currdate < rhdate && eqflag == 1){
			plot fMSnamefunc(i)				@YYEXIT1 lc rgb color(i),\
				fMSnamefunc(i)				@YCEXIT1 lc rgb color(i),\
				simpath.path(i,'')			@YYEXIT1SIM lc rgb color(i),\
				simpath.path(i,comflag)	@YYEXIT1SIMCOMP lc rgb color(i),\
				path_eq_vel(currdate)		@YYEXITEQ1VEL lc rgb color(i)
			}	
	}





		
	if (i < nexp) {
		set key off
		} else{
		set key inside right bottom
		}
}

unset title

#####################################################################################################
# Figure 2,1
#####################################################################################################
# @GFX1T
ifig = 1
jfig = 2

set origin xpos(ifig),ypost

@TYLABEL
@TXLABEL
@TYTICS
@TXTICS
@TSIZE

# add some data on the experiment to the plot 
y0 = -0.35 
diam = 10**(-1)  # in cm
sec = pi*(diam/2)**2
newy(n) = y0+n*y0/4.5

do for [i = 1:nexp] {
	set label at graph 0,newy((i-1)) fileid(i) font 'CMU Sans Serif,16' tc rgb color(i) noenhanced
}

set xrange [0:*]
set yrange [Tminh:Tmaxh]

do for [i=1:nexp] {
	set xrange [0:tmaxh]
	set yrange [Tminh:Tmaxh]
	set origin xpos(ifig),ypost
	@TSIZE
	
	if (i < nexp) {
		set tics format ''}
	{
		@TYTICS
		@TXTICS
	}	
	
	plot  	fLVnamefunc(i,'Tcol')		@TMID1 lc rgb color(i),\
			simpath.patht(i,'')			@TMID1SIM lc rgb color(i),\
			simpath.patht(i,comflag)	@TMID1SIMCOMP lc rgb color(i)
}

unset multiplot
set output

# open output
system 'start '.outfile