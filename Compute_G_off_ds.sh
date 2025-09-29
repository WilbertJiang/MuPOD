#!/bin/bash
for i in {1..11}; do
	echo $i
	for j in {1..11}; do
		python3 Compute_G_offdiag.py $i $j
	done
done
