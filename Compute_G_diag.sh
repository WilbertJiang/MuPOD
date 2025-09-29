#!/bin/bash
for i in {1..13}
do
	python3 Compute_G_diag_percell_dx.py $i & 
done
