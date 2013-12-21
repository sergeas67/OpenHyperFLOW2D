#!/bin/bash

HYPERFLOW2D=bin/OpenHyperFLOW2D-1.00
ProjectName=$1

MPI=`$HYPERFLOW2D | grep "parallel MPI"`

if [ "$MPI" != "" ]
then

if [ "$1" == "" ]
then
echo
echo "Usage: $0 {project name} [{num_nodes}]"
exit
fi

declare -i NHOSTS="$2"
if [ "${NHOSTS}" == "0" ]
then
NHOSTS=1
fi

declare -i CORESPERHOST=`cat /proc/cpuinfo | grep processor | wc -l`
declare -i NCORES=${NHOSTS}*${CORESPERHOST}
declare -i NUMNODES=${NHOSTS}+1
ulimit -s unlimited



`cat ./.mpi`/mpiexec -n ${NCORES} ${HYPERFLOW2D} ${ProjectName}.dat
#1>  ${ProjectName}.out  2> ${ProjectName}.err

else
if [ "$1" == "" ]
then
echo
echo "Usage: $0 {project name}"
exit
fi

${HYPERFLOW2D} ${ProjectName}.dat
#1>  ${ProjectName}.out  2> ${ProjectName}.err

fi