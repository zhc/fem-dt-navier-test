#!/bin/sh
#rm -rf result*
N=40
DT=1
T=30
CONV=d1
RE=100
mkdir png
for CONV in d1 d2 d3 n1 n2 n3 s1 s2 s3 s4 s5
#for CONV in d2 n2 s2 s3 s4 
#for CONV in d2 n2 s2
do
	python navier-dr.py ${N} ${DT} ${T} ${CONV} ${RE}
done 
