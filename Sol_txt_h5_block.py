# fenics code for the thermal simulation of chip for the publication  for frequency is 3.5 GHZ
from fenics import *
from dolfin import * 
from mshr import *
from numpy import loadtxt
from petsc4py import PETSc
import numpy as np
import sys
import time
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
#define the parameter
N_index = int(sys.argv[1])
N_STATE = 0                                                         # N_STATE is the controlling parameter   N_STATE = 0 indicates that we need to solve PDE and generate data; 	N_STATE = 1 indicate that we already have data just need to read it.
mesh_file_name = "./gmsh_mesh/Block" + str(N_index) + ".xdmf"
mesh_file = XDMFFile(comm, mesh_file_name)
mesh = Mesh()
mesh_file.read(mesh)
coor2 = mesh.coordinates()
####to generate the simulation mesh
lr = coor2[:,0].max()
ll = coor2[:,0].min()
wb = coor2[:,1].min()
wt = coor2[:,1].max()
hmax = coor2[:,2].max()
hmin = coor2[:,2].min()
T = 1200*3.125e-6                                                          # T is final time or total time (s)
V = FunctionSpace(mesh, 'P', 1)


num_steps = 1200                                                  # num_steps is the number of time steps
t = 0
dt = T/num_steps                                                # time step size 
Rb = 2.26                                                     # the thermal resistor of convectional face
#ch = 1/(Rb*Ach)                                                         # ch is convective coefficient
h_c = 2.40598e4
Ta = 0                                                         # Ta is initial /ambient temperature 
#Nu = 5                                                          # Nu is the number if functional unit
thick_actl =  0.0000557976                                        # lthick_actl is the total thickness of active layer (m)
data_filename = 'Block'+str(N_index)+'_sol_pred.txt'
rst = loadtxt(data_filename)
#import power trace file
#define initial value 
u0 = Constant(Ta)                                                                   # Ta is initial temperature 
u_n = interpolate(u0,V)
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
solution = []
for n in range(0,num_steps):
	solution.append('u1')
# Collect Robin integrals
#a = DS1*sc*u*v*dx + dt*kappa*dot(grad(u), grad(v))*dx + sum(integrals_R_a)
v2d = vertex_to_dof_map(V)
dd = 0
u = Function(V)
u = interpolate(u0,V)
for n in range(0,num_steps):# changed to save time
    print(dd)
    u = interpolate(u0,V)
    for i, kk in enumerate(coor2):
        u.vector()[v2d[i]]= rst[n][i]
    dd += 1
    solution_file_name = "solution_block"+str(N_index)+"_240_pred/file_" + str(n) + "h5"
    solution_file = HDF5File(mesh.mpi_comm(), solution_file_name, "w")
    solution_file.write(u,"solution")
    solution_file.close()
