#!/bin/bash
for i in {1..7..1}
do

	cp config_block$i.txt config_block.txt 
	./Prediction.sh $1 $i
done
