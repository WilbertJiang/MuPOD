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
T = 240*3.125e-6                                                          # T is final time or total time (s)
num_steps = 240                                                  # num_steps is the number of time steps
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
mesh = BoxMesh(Point(0.001220,0.007480,0), Point(0.008000,0.010930,0.000242),70,29,13)
V = FunctionSpace(mesh, 'P', 1)
Num_nodes = mesh.num_vertices()
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
# if we have the solutiuon, we just need to  Load solution# #######################read the solution from a solution file ######################
u = Function(V)
for n in range(0,Train_steps):
	solution_load_file_name = "../solution_block_2_240_pred/file_" + str(n) + "h5"
	solution_file = HDF5File(mesh.mpi_comm(), solution_load_file_name, "r")
	solution_file.read(u, "solution")
	solution[n] = interpolate(u0,V)
	solution[n].assign(u)
solution_file.close()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~generate podmode~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~generate podmode~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ################# generate the data matrix#####################
# ----write coorelation matrix into a text file-------
data_filename = './buidling_blk_2blocks/corelationmatrix.txt' 
mm = loadtxt(data_filename) 
co_matrix = PETScMatrix()   # co_matrix is correlation matrix
co_matrix.mat().setSizes([Train_steps,Train_steps])
co_matrix.mat().setType("dense")
co_matrix.mat().setUp()
for i in range(0,Train_steps):
	for j in range(0,Train_steps):
		if i <= j:
			co_matrix.mat().setValues(i,j,mm[i][j])
		else :
			co_matrix.mat().setValues(i,j,mm[j][i])
co_matrix.mat().assemble()
# ######################### solve eigenvalue and eigenvector############################
e_value_r = []
e_vec_r = []
# Create eigensolver
PETScOptions.set("st_ksp_type", "preonly")
PETScOptions.set("st_pc_type", "lu")
PETScOptions.set(" -pc_factor_mat_solver_type", "mumps")
eigensolver = SLEPcEigenSolver(co_matrix)
# Configure solver to get the largest eigenvalues first:
eigensolver.parameters["spectrum"] = "largest real"
# Compute all eigenvalues of A x = \lambda x
print ('Computing eigenvalues. This can take a minute')
print(' solve %d largest eigenpairs' %(Train_steps))
eigensolver.solve()
header = './buidling_blk_2blocks/eigenvalue_r'                               
eigenvalue_file_name = header + '.txt'
eigenvalue_file = open(eigenvalue_file_name,'w')
for i in range (0,Train_steps):
	# Extract largest (first) eigenpair
	r, c, rx , cx = eigensolver.get_eigenpair(i)
	e_value_r.append(r)
	e_vec_r.append(rx)
	eigenvalue_file.write('%.16g\n' % (r))
	#print(' print the eigenvalue \n')
	#print('%8g\n' %(r))
eigenvalue_file.close()
print('eigenvalue is done \n')
# ###########################################################################
# ###########################################################################
# ######################### generate the pod mode ############################
# ###########################################################################
# ###########################################################################
podmode = []

#podmode_vector_file = open('podmode_fenics_core1_test_square_lp_redo.txt','w')      # save the podmode into a file 
for i in range (0,Train_steps):
#for i in range (0,Train_steps):
    a_init = Constant(0)
    a = interpolate(a_init,V)
    for j in range (0, Train_steps):	
        a.vector().axpy(e_vec_r[i][j], solution[j].vector())
    a.vector()[:] =  a.vector()/(e_value_r[i]*Train_steps)
	# **********normalize the podmode ***************
    normalize_value = assemble(dot(a,a)*dx)
    a.vector()[:] =  a.vector()/sqrt(normalize_value)
    #print(assemble(dot(a,a)*dx))
    podmode.append(a)
    pod_filename = './buidling_blk_2blocks/podmode_fenics.txt'
    podmode_vector_file = open(pod_filename,'a')
    for n_count in range(0,Num_nodes):
        podmode_vector_file.write('%.16g\t' % (a.vector()[v2d[n_count]]))
    podmode_vector_file.write('\n')
    podmode_vector_file.close( )
    mode_file_name = "./buidling_blk_2blocks/podmode/mode_" + str(i) + "h5"
    mode_file = HDF5File(mesh.mpi_comm(), mode_file_name, "w")
    mode_file.write(a,"solution")
    mode_file.close()


#print(assemble(dot(a,a)*dx))
#######check
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ save podmode into file~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ save podmode into file~~~~~~~~~~~~~~~~~~~~~~~

