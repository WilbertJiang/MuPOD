# fenics code for the thermal simulation of chip for the publication  for frequency is 3.5 GHZ
import sys
import meshio
from fenics import *
from dolfin import * 
from mshr import *
from numpy import loadtxt
from petsc4py import PETSc
import numpy as np
import time
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
#define the parameter													# w is the width of domain m
h = 0.0002417896															# h is the thickness of domain m												# hs is number of grid along with h direction
T = 240*3.125e-6                                                          # T is final time or total time (s)
num_steps = 240                                                  # num_steps is the number of time steps
t = 0
Train_steps = num_steps
dt = T/num_steps                                                # time step size 
Rb = 2.26                                                        # the thermal resistor of convectional face
#ch = 1/(Rb*Ach)                                                         # ch is convective coefficient
h_c = 2.40598e4
Ta = 0     
config_block = loadtxt('config_block.txt',dtype="int")
# Ta is initial /ambient temperature 
N_mode = config_block[0][1]                                                      # Nu is the number if functional unit
#N_mode = config_block[N_index-1][1]                                                       # Nu is the number if functional unit
flp = loadtxt("../Floorplan_AMD_multiblock_cutting.txt")
                                               #chip_area is the area of chip
# ########################################################################################################################################
#                              import geometry and mesh from Gmsh                                                                        #
# ########################################################################################################################################
# import mesh file ---------multiblock.msh
#mesh_file = XDMFFile(comm, "mesh1.xdmf")
#mesh = Mesh()
#mesh_file.read(mesh)
#M = Box(Point(0,0,0), Point(l,w,h))
mesh = BoxMesh(Point(0.001220,0.007480,0), Point(0.008000,0.010930,0.000242),70,29,13)
V = FunctionSpace(mesh, 'P', 1)
Num_nodes = mesh.num_vertices()
thick_actl =  0.0000557976
#number_mode = loadtxt('config_block.txt')
#define thermal conductivity 
tol = 1E-14
k_0 = 100														# silicon conductivity      (w/(m.k))
k_1 = 100                                                      #oxide silicon thermal conductivity
kappa = Expression('x[2] <= 0.00045 + tol ? k_0 : k_1', degree=0,tol=tol, k_0=k_0, k_1=k_1) #define subdomain
#define density 
D0 = 2.33e3                                                         # silicon density    (kg/m^3)
D1 = 2.33e3                                                        # oxide silicon density
DS1 = Expression('x[2] <= 0.00045 + tol ? D0 : D1', degree=0,tol=tol, D0=D0, D1=D1) #define subdomain
#define specific heat
c0 = 751.1                                                         # silicon specific heat   (J/(kg.k))
c1 = 751.1                                                         # oxide silicon specific heat
sc = Expression('x[2] <= 0.00045 + tol ? c0 : c1', degree=0,tol=tol, c0=c0, c1=c1) #define subdomain
#define power source term 

T_integral = 0
#define initial value 
u0 = Constant(Ta)                                                                   # Ta is initial temperature 
u_n = interpolate(u0,V)
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
solution = []
for n in range(0,Train_steps):
	solution.append('u1')
# Collect Neumann integrals
#a = DS1*sc*u*v*dx + dt*kappa*dot(grad(u), grad(v))*dx + sum(integrals_R_a)
coor = mesh.coordinates()
v2d = vertex_to_dof_map(V)
h = coor[:,2].max()
CU = loadtxt('pod_result/CU.txt')
#print(len(CU[0]))
# if we have the solutiuon, we just need to  Load solution# #######################read the solution from a solution file ######################
u = Function(V)
f = interpolate(u0,V)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~generate podmode~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ###########################################################################
"""
for i in range(0,Num_nodes):
    if abs(coor[i][2] - coor[:,2].max())< tol and abs(coor[i][0] - 0.00361) < 0.00012:
        print(i)
"""
podmode = []
for n in range(0,10):
    podmode.append('u1')
for n in range(0,8):
    podmode_load_file_name = "/home/jiangl2/multi_block_AMD_240_cutting_smalldiff/Building_block/buidling_blk_2blocks/podmode/mode_" + str(n) + "h5"
    podmode_file = HDF5File(mesh.mpi_comm(), podmode_load_file_name, "r")
    podmode_file.read(u, "solution")
    podmode[n] = interpolate(u0,V)
    podmode[n].assign(u)
    podmode_file.close()
################compute the index of coefficient
id_total = 0
#index = 3460
#######check
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ save podmode into file~~~~~~~~~~~~~~~~~~~~~~~
N_index = 9

#import matplotlib.pyplot as plt
for p_i in range(0,8):
    N_index2 = p_i
    tol = 1e-14 # avoid hitting points outside the domain
    y = np.linspace(0.007480+ tol, 0.010930- tol, 100)
#points = [(0.00239 + 0.00122, y_,h) for y_ in y]  # 2D points
    points = [(0.00239 + 0.00122, y_,h) for y_ in y]  # 2D points

    p_line = np.array([podmode[p_i](point) for point in points])

    header_data = 'podmode_plot/y_block'+str(N_index)
    data_file_name = header_data + '.txt'
    data_file = open(data_file_name,'w')
    for i in range(0,len(y)):
        data_file.write('%.16g\n' % (y[i]))
    data_file.close()

    header_data = 'podmode_plot/T_Y_block'+str(N_index)+'_'+str(N_index2)+'mode'
    data_file_name = header_data + '.txt'
    data_file = open(data_file_name,'w')
    for i in range(0,len(y)):
        data_file.write('%.16g\n' % (p_line[i]))
    data_file.close()
    
    tol = 1e-14
    y = np.linspace(0.001220 + tol, 0.008000- tol, 100)
    points = [(y_,0.009205,h) for y_ in y]
    p_line = np.array([podmode[p_i](point) for point in points])
    header_data = 'podmode_plot/x_block'+str(N_index)
    data_file_name = header_data + '.txt'
    data_file = open(data_file_name,'w')
    for i in range(0,len(y)):
        data_file.write('%.16g\n' % (y[i]))
    data_file.close()
    header_data = 'podmode_plot/T_X_block'+str(N_index)+'_'+str(N_index2)+'mode'
    data_file_name = header_data + '.txt'
    data_file = open(data_file_name,'w')
    for i in range(0,len(y)):
        data_file.write('%.16g\n' % (p_line[i]))
    data_file.close()
#######generate the pvd file

