#!/bin/bash
for i in {1..13}
do
	#echo $i
	for j in {1..13}; do
		echo $i $j
		python3 Compute_G_offdiag_percell.py $i $j& 
	done
done
