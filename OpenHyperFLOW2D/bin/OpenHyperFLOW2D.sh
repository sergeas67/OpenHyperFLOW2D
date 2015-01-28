#!/bin/bash

HYPERFLOW2D=bin/OpenHyperFLOW2D-1.02
ProjectName=$1
MPILIB=`cat ./.mpilib`
MPI=`LD_LIBRARY_PATH=${MPILIB} $HYPERFLOW2D | grep "parallel MPI"`

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

declare -i CORESPERHOST=`cat /proc/cpuinfo | grep processor | wc -l`/2
declare -i NCORES=${NHOSTS}*${CORESPERHOST}
declare -i NUMNODES=${NHOSTS}+1
ulimit -s unlimited


# MPI 1
#declare -i pnum=0
#rm -f ./.hosts
#while [ $pnum -lt $CORESPERHOST ]
#do
#echo localhost >> ./.hosts
#pnum=$pnum+1
#done
#LD_LIBRARY_PATH=${MPILIB} `cat ./.mpi`/mpirun --bind-to-core  -np ${NCORES} -hostfile .hosts ${HYPERFLOW2D} ${ProjectName}.dat
#-
# MPI 2
LD_LIBRARY_PATH=${MPILIB}  `cat ./.mpi`/mpiexec -n ${NCORES}  ${HYPERFLOW2D} ${ProjectName}.dat
# 1>  ${ProjectName}.out  2> ${ProjectName}.err
#  taskset -c 0,2,4,6  

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