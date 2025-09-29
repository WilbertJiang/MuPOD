#!/bin/bash
for i in {1..13} 
do
	mpirun -n 20 python3 ComputingAM_flp7.py $i
	python3 save_podmode.py $i
done
