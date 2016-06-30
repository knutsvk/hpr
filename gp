#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.4 patchlevel 3
#    	last modified March 2011
#    	System: Linux 3.13.0-83-generic
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2010
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help seeking-assistance"
#    	immediate help:   type "help"
#    	plot window:      hit 'h'
# set terminal wxt 0
# set output
unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border 31 front linetype -1 linewidth 1.000
set xdata
set ydata
set zdata
set x2data
set y2data
set timefmt x "%d/%m/%y,%H:%M"
set timefmt y "%d/%m/%y,%H:%M"
set timefmt z "%d/%m/%y,%H:%M"
set timefmt x2 "%d/%m/%y,%H:%M"
set timefmt y2 "%d/%m/%y,%H:%M"
set timefmt cb "%d/%m/%y,%H:%M"
set boxwidth
set style fill  empty border
set style rectangle back fc lt -3 fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0, 0 
set dummy x,y
set format x "% g"
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set format cb "% g"
set angles radians
unset grid
set key title ""
set key inside right top vertical Right noreverse enhanced autotitles nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
unset label
unset arrow
set style increment default
unset style line
set style line 1  linetype 1 linecolor rgb "#352a87"  linewidth 1.000 pointtype 1 pointsize default pointinterval 0
set style line 2  linetype 1 linecolor rgb "#0f5cdd"  linewidth 1.000 pointtype 2 pointsize default pointinterval 0
set style line 3  linetype 1 linecolor rgb "#1481d6"  linewidth 1.000 pointtype 3 pointsize default pointinterval 0
set style line 4  linetype 1 linecolor rgb "#06a4ca"  linewidth 1.000 pointtype 4 pointsize default pointinterval 0
set style line 5  linetype 1 linecolor rgb "#2eb7a4"  linewidth 1.000 pointtype 5 pointsize default pointinterval 0
set style line 6  linetype 1 linecolor rgb "#87bf77"  linewidth 1.000 pointtype 6 pointsize default pointinterval 0
set style line 7  linetype 1 linecolor rgb "#d1bb59"  linewidth 1.000 pointtype 7 pointsize default pointinterval 0
set style line 8  linetype 1 linecolor rgb "#fec832"  linewidth 1.000 pointtype 8 pointsize default pointinterval 0
set style line 9  linetype 1 linecolor rgb "#f9fb0e"  linewidth 1.000 pointtype 9 pointsize default pointinterval 0
set style line 11  linetype 1 linecolor rgb "#0072bd"  linewidth 1.000 pointtype 11 pointsize default pointinterval 0
set style line 12  linetype 1 linecolor rgb "#d95319"  linewidth 1.000 pointtype 12 pointsize default pointinterval 0
set style line 13  linetype 1 linecolor rgb "#edb120"  linewidth 1.000 pointtype 13 pointsize default pointinterval 0
set style line 14  linetype 1 linecolor rgb "#7e2f8e"  linewidth 1.000 pointtype 14 pointsize default pointinterval 0
set style line 15  linetype 1 linecolor rgb "#77ac30"  linewidth 1.000 pointtype 15 pointsize default pointinterval 0
set style line 16  linetype 1 linecolor rgb "#4dbeee"  linewidth 1.000 pointtype 16 pointsize default pointinterval 0
set style line 17  linetype 1 linecolor rgb "#a2142f"  linewidth 1.000 pointtype 17 pointsize default pointinterval 0
unset style arrow
set style histogram clustered gap 2 title  offset character 0, 0, 0
unset logscale
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view map
set samples 100, 100
set isosamples 10, 10
unset surface
unset contour
set clabel '%8.3g'
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 1 1,1
set origin 0,0
set style data points
set style function lines
set xzeroaxis linetype -2 linewidth 1.000
set yzeroaxis linetype -2 linewidth 1.000
set zzeroaxis linetype -2 linewidth 1.000
set x2zeroaxis linetype -2 linewidth 1.000
set y2zeroaxis linetype -2 linewidth 1.000
set ticslevel 0.5
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set xtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0
set xtics autofreq  norangelimit
set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0
set ytics autofreq  norangelimit
set ztics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0
set ztics autofreq  norangelimit
set nox2tics
set noy2tics
set cbtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0
set cbtics autofreq  norangelimit
set title "" 
set title  offset character 0, 0, 0 font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  offset character 0, 0, 0 font "" norotate
set rrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set trange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set urange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set vrange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set xlabel "" 
set xlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set x2label "" 
set x2label  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set xrange [ * : * ] noreverse nowriteback  # (currently [0.00000:1.00000] )
set x2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set ylabel "" 
set ylabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set y2label "" 
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set yrange [ * : * ] noreverse nowriteback  # (currently [1.00000:0.00000] )
set y2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set zlabel "" 
set zlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse nowriteback  # (currently [-1.00000:1.00000] )
set cblabel "" 
set cblabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set cbrange [ * : * ] noreverse nowriteback  # (currently [-1.00000:1.00000] )
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "en_GB.UTF-8"
set pm3d implicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette defined ( 0 0.2078 0.1647 0.5294, 0.125 0.01176 0.3882 0.8824, 0.25 0.07843 0.5216 0.8314,\
     0.375 0.02353 0.6549 0.7765, 0.5 0.2196 0.7255 0.6196, 0.625 0.5725 0.749 0.451, 0.75 0.851 0.7294 0.3373,\
     0.875 0.9882 0.8078 0.1804, 1 0.9765 0.9843 0.0549 )
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set loadpath 
set fontpath 
set fit noerrorvariables
#    EOF
