#!/bin/bash
HYPERFLOW2D=bin/OpenHyperFLOW2D-CUDA-1.02
ProjectName=$1
if [ "$1" == "" ]
then
echo
echo "Usage: $0 {project name}"
exit
fi

#optirun --no-xorg

${HYPERFLOW2D} ${ProjectName}.dat

#1>  ${ProjectName}.out  2> ${ProjectName}.err
