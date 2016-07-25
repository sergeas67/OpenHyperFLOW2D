#!/bin/bash
TERM_TYPE="x11"
EXT="png"

PLT_FILE1=$1
PLT_FILE2=$2
PLT_FILE3=$3
PLT_FILE4=$4

head -n1  $PLT_FILE1 | sed s/"#VARIABLES = "/""/g | tr -s "," '\012' > var.txt
declare -i NUM_PARAM=`cat var.txt | wc -l `+1

declare -i INDEX=2

while  [ $INDEX -lt $NUM_PARAM ] 
do

TITLE=`head -n $INDEX var.txt | tail -n1`

IS_2D=`echo $TITLE | grep "(N)"`
#IS_2D=`echo $TITLE | grep "(Dp)"`
IS_3D=`echo $TITLE | grep "(x;y)"`
TITLE=`echo $TITLE | tr -s ";" ","`

cat /dev/null >  $TITLE.gplt

if [ "$TERM_TYPE" == "latex" ]
then
cat /dev/null >  _$TITLE.$EXT
cat /dev/null >  $TITLE.$EXT
fi

if [ $IS_3D ]
then
 PLOT_CMD="splot '$PLT_FILE1' u 1:2:$INDEX  t '$TITLE'"
else
 
 if [ "$PLT_FILE1" ]
 then
  PLOT_CMD="plot '$PLT_FILE1' u 1:(\$$INDEX>0 ? \$$INDEX : 1/0)  t '$TITLE'"

  if [ "$PLOT_ALL" == "" ] 
  then
  PLOT_ALL="plot '$PLT_FILE1' u 1:(\$$INDEX>0 ? \$$INDEX : 1/0)  t '$TITLE'"
  else
  PLOT_ALL+=", '' u 1:(\$$INDEX>0 ? \$$INDEX : 1/0)  t '$TITLE'"
  fi
  
 fi
 
 if [ "$PLT_FILE2" ]
 then
 PLOT_CMD="$PLOT_CMD ,'$PLT_FILE2' u 1:(\$$INDEX>0 ? \$$INDEX : 1/0)  t '$TITLE'"
 fi

 if [ "$PLT_FILE3" ]
 then
 PLOT_CMD="$PLOT_CMD ,'$PLT_FILE3' u 1:(\$$INDEX>0 ? \$$INDEX : 1/0)  t '$TITLE'"
 fi

 if [ "$PLT_FILE4" ]
 then
 PLOT_CMD="$PLOT_CMD ,'$PLT_FILE4' u 1:(\$$INDEX>0 ? \$$INDEX : 1/0)  t '$TITLE'"
 fi

fi


if [ "$TERM_TYPE" == "latex" ]
then
cat >> $TITLE.$EXT << _TEX
\documentclass[a4paper,10pt]{article}
\usepackage[T1]{fontenc}
\usepackage[koi8-r]{inputenc}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{latexsym}
\pdfcompresslevel=9
\begin{document}                  
\input{_$TITLE.$EXT}
\end{document}
_TEX
fi

if [ "$TERM_TYPE" != "x11" ]
then
SET_TERM="set terminal $TERM_TYPE"
else
SET_TERM="set terminal $TERM_TYPE $INDEX"
fi

if [ "$TERM_TYPE" != "x11" ]
then
SET_OUTPUT="set output \"_$TITLE.$EXT\""
else
PAUSE="pause 10"
fi

IS_RMS=`echo $TITLE | grep RMS`

if [ "$IS_RMS" != "" ]
then
LOGSCALE="set logscale y 10"
else
LOGSCALE="#"
fi


cat >> $TITLE.gplt << CMD
$SET_TERM

set view 0, 0, 1
set nosurface
set contour
set data style lines
#histeps
#lines
#set yrange [-200:]
set cntrparam levels auto 40
set grid
set title "$TITLE"
$LOGSCALE
set logscale y 10
show title
$SET_OUTPUT
$PLOT_CMD
$PAUSE
CMD

#gnuplot  $TITLE.gplt

if [ "$TERM_TYPE" == "latex" ]
then
pdflatex $TITLE.$EXT
rm -f    $TITLE.$EXT _$TITLE.$EXT $TITLE.aux $TITLE.log
fi
rm -f $TITLE.gplt
INDEX=$INDEX+1;
done
rm -f RMS_All.gplt
cat >> RMS_All.gplt << CMD_ALL
$SET_TERM

set view 0, 0, 1
set nosurface
set contour
set style data lines
set cntrparam levels auto 40
set grid
set title "All Residuals RMS"
#set yrange [1.e-8:]
$LOGSCALE
set logscale y 10
show title
$SET_OUTPUT
$PLOT_ALL
$PAUSE
#pause -1
CMD_ALL

gnuplot  RMS_All.gplt
#rm -f RMS_All.gplt
