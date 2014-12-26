#!/bin/bash
#LD_LIBRARY_PATH=${MPILIB}
HYPERFLOW2D=bin/OpenHyperFLOW2D-1.02
ProjectName=$1
MPILIB=`cat ./.mpilib`
MPI=`$HYPERFLOW2D | grep parallel | grep MPI`

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


# MPI 1
#declare -i pnum=0
#rm -f ./.hosts
#while [ $pnum -lt $CORESPERHOST ]
#do
#echo localhost >> ./.hosts
#pnum=$pnum+1
#done
#LD_LIBRARY_PATH=${MPILIB} `cat ./.mpi`/mpirun -np ${NCORES} -hostfile .hosts ${HYPERFLOW2D} ${ProjectName}.dat

# MPI 2
LD_LIBRARY_PATH=${MPILIB}  `cat ./.mpi`/mpiexec -n ${NCORES}  ${HYPERFLOW2D} ${ProjectName}.dat 
#1>  ${ProjectName}.out  2> ${ProjectName}.err

else
if [ "$1" == "" ]
then
echo
echo "Usage: $0 {project name}"
exit
fi

#optirun --no-xorg 
${HYPERFLOW2D} ${ProjectName}.dat
#1>  ${ProjectName}.out  2> ${ProjectName}.err

fi