#!/bin/bash

FLAGS=`cat /proc/cpuinfo | grep flags | sort -u | cut -d : -f 2`

SIMD=None

for FLAG in $FLAGS
do

if [ "$FLAG" == "sse" ]
then
SIMD=SSE
fi

if [ "$FLAG" == "sse2" ]
then
SIMD=SSE2
fi

if [ "$FLAG" == "ssse3" ]
then
SIMD=SSE3
fi

if [ "$FLAG" == "sse4_1" ]
then
SIMD=SSE41
fi

if [ "$FLAG" == "sse4_2" ]
then
SIMD=SSE42
fi

if [ "$FLAG" == "avx" ]
then
SIMD=AVX
fi


done
echo -e $SIMD