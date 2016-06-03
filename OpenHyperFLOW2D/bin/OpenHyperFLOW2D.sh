#!/bin/bash

HYPERFLOW2D=bin/OpenHyperFLOW2D-1.03
ProjectName=$1

MPILIB=`cat ./.mpilib`

declare -i NHOSTS="$2"
if [ "${NHOSTS}" == "0" ]
then
NHOSTS=1
fi

declare -i CORESPERHOST=`cat /proc/cpuinfo | grep processor | wc -l`
declare -i NCORES=${NHOSTS}*${CORESPERHOST}

#echo ${NCORES}

MPI=`LD_LIBRARY_PATH=${MPILIB} $HYPERFLOW2D | grep "parallel MPI"`

if [ "$MPI" != "" ]
then

if [ "$1" == "" ]
then
echo
echo "Usage: $0 {project name} [{num_nodes}]"
exit
fi


declare -i NUMNODES=${NHOSTS}+1
ulimit -s unlimited


# MPI 2
LD_LIBRARY_PATH=${MPILIB}  `cat ./.mpi`/mpiexec -np ${NCORES}  ${HYPERFLOW2D} ${ProjectName}.dat
# 1>  ${ProjectName}.out  2> ${ProjectName}.err

else
if [ "$1" == "" ]
then
echo
echo "Usage: $0 {project name}"
exit
fi

OMP_NUM_THREADS=${NCORES}   ${HYPERFLOW2D} ${ProjectName}.dat
#1>  ${ProjectName}.out  2> ${ProjectName}.err

fi