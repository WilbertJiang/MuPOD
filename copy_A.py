# fenics code for the thermal simulation of chip for the publication  for frequency is 3.5 GHZ
import meshio
from fenics import *
from dolfin import * 
from mshr import *
from numpy import loadtxt
from petsc4py import PETSc
import sys
import numpy as np
import time
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
#define the parameter													# w is the width of domain m
h = 0.0002417896															# h is the thickness of domain m												# hs is number of grid along with h direction
T = 1200*4.347826086956521e-6                                                          # T is final time or total time (s)
num_steps = 1200                                                  # num_steps is the number of time steps
t = 0
Train_steps = num_steps
dt = T/num_steps                                                # time step size 
Rb = 2.26                                                        # the thermal resistor of convectional face
#ch = 1/(Rb*Ach)                                                         # ch is convective coefficient
h_c = 2.40598e4
Ta = 0                                                         # Ta is initial /ambient temperature 
N_mode = 30                                                       # Nu is the number if functional unit
                                               #chip_area is the area of chip
# ########################################################################################################################################
#                              import geometry and mesh from Gmsh                                                                        #
# ########################################################################################################################################
# import mesh file ---------multiblock.msh
#mesh_file = XDMFFile(comm, "mesh1.xdmf")
#mesh = Mesh()
#mesh_file.read(mesh)
#M = Box(Point(0,0,0), Point(l,w,h))
N_index = int(sys.argv[1])
# ----write coorelation matrix into a text file-------
data_filename = '/home/jiangl2/multi_block_AMD_240_cutting_fixed_small_1800_realpower/Building_block/buidling_blk_'+ str(N_index)+'/corelationmatrix.txt' 
mm = loadtxt(data_filename) 

header = './buidling_blk_'+ str(N_index)+'/corelationmatrix'
C_matrix_file_name = header + '.txt'
C_matrix_file = open(C_matrix_file_name,'w')
for i in range (0,num_steps):
        for j in range (0,num_steps):
                C_matrix_file.write('%.16g\t' % (1800*mm[i][j]/1200))
        C_matrix_file.write('\n')
C_matrix_file.close()
