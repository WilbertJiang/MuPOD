#!/bin/bash
./Construct_C
./Construct_G $1
./Construct_P
./ODE_solver $[13*$2]
for i in {1..13..1}
do
	python3 Prediction_MLB_CPU.py $i $i
done
