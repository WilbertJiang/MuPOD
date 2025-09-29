# fenics code for the thermal simulation of chip for the publication  for frequency is 3.5 GHZ
import meshio
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
mesh_file_name = "/home/jiangl2/multi_block_AMD_1200/gmsh_mesh/Block" + str(N_index) + ".xdmf"
mesh_file = XDMFFile(comm, mesh_file_name)
mesh = Mesh()
mesh_file.read(mesh)

N_block = N_index

Num_nodes = mesh.num_vertices()
V = FunctionSpace(mesh, 'P', 1)
coor = mesh.coordinates()
v2d = vertex_to_dof_map(V)

lr = coor[:,0].max()
ll = coor[:,0].min()
wb = coor[:,1].min()
wt = coor[:,1].max()
hmax = coor[:,2].max()
hmin = coor[:,2].min()
T = 1200*3.125e-6                                                             # T is final time or total time (s)
num_steps = 1200            
Train_steps = 1200# num_steps is the number of time steps
t = 0
dt = T/num_steps                                                # time step size 
Rb = 2.26                                                     # the thermal resistor of convectional face
#ch = 1/(Rb*Ach)                                                         # ch is convective coefficient
h_c = 2.40598e4
Ta = 0                                                         # Ta is initial /ambient temperature 
N_mode = 30
Nu = 11                                                        # Nu is the number if functional unit
thick_actl = 0.0000557976                                         # lthick_actl is the total thickness of active layer (m)
# ########################################################################################################################################
#                              import geometry and mesh from Gmsh                                                                        #
# ########################################################################################################################################
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

# Define boundary subdomains
class BoundaryX0(SubDomain):
      def inside(self, x, on_boundary):
          return on_boundary and near(x[0], ll, tol)

class BoundaryX1(SubDomain):
    def inside(self, x, on_boundary):
          return on_boundary and near(x[0], lr, tol)

class BoundaryY0(SubDomain):
      def inside(self, x, on_boundary):
          return on_boundary and near(x[1], wb, tol)

class BoundaryY1(SubDomain):
      def inside(self, x, on_boundary):
          return on_boundary and near(x[1], wt, tol)
class BoundaryZ0(SubDomain):
      def inside(self, x, on_boundary):
          return on_boundary and near(x[2], hmax, tol)

class BoundaryZ1(SubDomain):
      def inside(self, x, on_boundary):
          return on_boundary and near(x[2], hmin, tol)
boundary_markers = MeshFunction('size_t', mesh,2)
boundary_markers.set_all(9999)
bx0 = BoundaryX0()
bx1 = BoundaryX1()
by0 = BoundaryY0()
by1 = BoundaryY1()
bz0 = BoundaryZ0()
bz1 = BoundaryZ1()
bx0.mark(boundary_markers, 0)
bx1.mark(boundary_markers, 1)
by0.mark(boundary_markers, 2)
by1.mark(boundary_markers, 3)
bz0.mark(boundary_markers, 4)
bz1.mark(boundary_markers, 5)
# Redefine boundary integration measure
ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
T_integral = 0
# Define boundary conditions
boundary_conditions = {0: {'Neumann': 0},   # x = 0 adiabatic boundary condition
                       1: {'Neumann': 0},   # x = 1
                       2: {'Neumann': 0},   # y = 0
                       3: {'Neumann': 0},    # y = 1
		       4: {'Neumann': 0}, # z = 0
		       5: {'Robin': (h_c, Ta)}}      # z = 1       r is dt*h , s = Ta reference temperature
					   #5:{'Robin': (dt*k_0, Ta)}}      # z = 1       r is dt*h , s = Ta reference temperature

# Collect Dirichlet conditions
bcs = []
for i in boundary_conditions:
	if 'Dirichlet' in boundary_conditions[i]:
		bc = DirichletBC(V, boundary_conditions[i]['Dirichlet'], boundary_markers, i)
		bcs.append(bc) 
#define initial value 
u0 = Constant(Ta)                                                                   # Ta is initial temperature 
u_n = interpolate(u0,V)
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
podmode = []

for n in range(0,Train_steps):
	podmode.append('u1')
solution = []
for n in range(0,Train_steps):
    solution.append('u1')
# Collect Neumann integrals
integrals_N = []
for i in boundary_conditions:
    if 'Neumann' in boundary_conditions[i]:
        if boundary_conditions[i]['Neumann'] != 0:
            g = boundary_conditions[i]['Neumann']
            integrals_N.append(g*v*ds(i))	
# Collect Robin integrals
integrals_R_a = []
integrals_R_L = []
for i in boundary_conditions:
	if 'Robin' in boundary_conditions[i]:
		r, s = boundary_conditions[i]['Robin']
		integrals_R_a.append(r*u*v*ds(i))
		integrals_R_L.append(r*s*v*ds(i))			

# Sum integrals to define variational problem
#a = DS1*sc*u*v*dx + dt*kappa*dot(grad(u), grad(v))*dx + sum(integrals_R_a)
coor = mesh.coordinates()
v2d = vertex_to_dof_map(V)
# if we have the solutiuon, we just need to  Load solution# #######################read the solution from a solution file ######################
u = Function(V)
for n in range(0,Train_steps):
    podmode_load_file_name = "/home/jiangl2/multi_block_AMD_1200/Building_block/buidling_blk_"+str(N_block)+"/podmode/mode_" + str(n) + "h5"
    #podmode_load_file_name = "./podmode_flp7/mode_" + str(n) + "h5"
    podmode_file = HDF5File(mesh.mpi_comm(), podmode_load_file_name, "r")
    podmode_file.read(u, "solution")
    #print(assemble(u*dx))
    podmode[n] = interpolate(u0,V)
    podmode[n].assign(u)
    podmode_file.close()
#print(assemble(podmode[0]*ds(4)))
for n in range(0,Train_steps):
    #solution_load_file_name = "../solution_flp7_240/file_" + str(n) + "h5"
    solution_load_file_name = "/home/jiangl2/multi_block_AMD_1200/solution_block"+str(N_block)+"_240_pred/file_" + str(n) + "h5"
    solution_file = HDF5File(mesh.mpi_comm(), solution_load_file_name, "r")
    solution_file.read(u, "solution")
    solution[n] = interpolate(u0,V)
    solution[n].assign(u)
    solution_file.close()
###########################compute the diagonal-Boundary term of G########################
n = FacetNormal(mesh)
header = "temp_gradient_lib_block"+str(N_block)
temp_gradient_file_name = header +'.txt'
temp_gradient_file = open(temp_gradient_file_name,'w')
temp_gradient = np.zeros((N_mode, num_steps))
for i in range(0, N_mode):
    for j in range(0, num_steps):
        temp_gradient[i][j]=assemble(kappa*dot(podmode[i],inner(grad(solution[j]),n))*ds(4))+assemble(kappa*dot(podmode[i],inner(grad(solution[j]),n))*ds(5))
        temp_gradient_file.write('%.16g\t'%(temp_gradient[i][j]))
    temp_gradient_file.write('\n')
temp_gradient_file.close()

u_n = interpolate(u0,V)
for i in range(0,Num_nodes):
    if coor[i][2] > hmax - thick_actl -tol:
        u_n.vector()[v2d[i]] =1

P_matrix = np.zeros((N_mode, 1))
header = "P_lib_block"+str(N_block)
P_matrix_file_name = header + '.txt'
P_matrix_file = open(P_matrix_file_name,'w')
for i in range (0,N_mode):
    P_matrix[i][0] = assemble(dot(podmode[i],u_n)*dx)
    P_matrix_file.write('%.16g\n' % (P_matrix[i][0]))
P_matrix_file.close()



