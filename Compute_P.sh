#!/bin/bash
for i in {1..13}; do
	cd /home/jiangl2/multi_block_AMD_240_cutting_fixed_small_1200_realpower_hotspot/Building_block/buidling_blk_$i
	python3 Compute_P_$i.py $i & 	
done
