#!/bin/bash
for i in {1..13}
do
	echo $i
	python3 Sol_simu_block.py $i
	python3 Sol_txt_h5_block.py $i
done	
