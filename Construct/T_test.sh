#!/bin/bash
cp config_block$2.txt config_block.txt 
./Construct_C
./Construct_G $1
./Construct_P
./ODE_solver $[13*$2]
python3 Prediction_MLB_CPU.py 9 9
