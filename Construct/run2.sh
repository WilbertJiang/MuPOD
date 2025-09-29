#!/bin/bash
./Construct_C

./Construct_G 40
./Construct_P
./ODE_solver $1
python3  Prediction_MLB_CPU.py 7 1
python3  Prediction_MLB_CPU.py 9 2
