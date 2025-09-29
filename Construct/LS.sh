#!/bin/bash
BASE=0.001
for counter1 in {1..5..1}
do
	for counter2 in {1..9..1}
	do

		a=$[$counter2*$[10**$counter1]]
		b=$(echo "$BASE*$a"|bc)
		echo $b
		./Prediction.sh $b $1
		python3 ./pod_result/Cal_error_num.py $1  
		rm ./pod_result/LS_error_each_block_$1.txt
		rm ./pod_result/T*.txt
		rm ./pod_result/x*.txt
		rm ./pod_result/y*.txt




	done
done
