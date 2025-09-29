#!/bin/bash
echo "FEM"
mpirun -n 20 python3 FEM_simulation_flp79_usingfenicsmesh.py
echo "write data"
python3 Sol_simu_block.py
echo "write block7"
python3  Sol_simu_block7.py
echo "generate h5"
python3 Sol_txt_h5_block.py
echo "generate block 5"
python3  Sol_txt_h5_block7.py

echo "go to block9"
cd ./buidling_blk_9
echo "Computing A"
mpirun -n 10 python3 ComputingAM_flp7.py
echo "save mode"
python3 save_podmode.py
echo "go to block 7"
cd ../buidling_blk_7
echo "Computing A"
mpirun -n 10 python3 ComputingAM_flp7.py
python3 save_podmode.py
echo "G_diag"
python3 Compute_G_diag_percell_ds.py
echo "G_offdiag"
python3 Compute_G_offdiag_percell_ds.py
echo "go to block 9"
cd ../buidling_blk_9
echo "G_diag"
python3 Compute_G_diag_percell.py
echo "G_offdiag"
python3  Compute_G_offdiag_percell_ds.py 
