#!/bin/bash
#echo Wait...15 sec
#sleep 15
if [ "$1" == "" ]
then
echo "Usage: $0 RMS_data_file"
else
rm -f RMS_All.gplt
./viewplt.sh $1

while true
do
gnuplot RMS_All.gplt
done
fi