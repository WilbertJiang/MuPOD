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
parms = parameters["krylov_solver"]
parms["relative_tolerance"]=1e-23
parms["absolute_tolerance"]=1e-25

mesh_file_name = "./gmsh_mesh/Block" + str(N_index) + ".xdmf"
mesh_file = XDMFFile(comm, mesh_file_name)
mesh = Mesh()
mesh_file.read(mesh)
coor2 = mesh.coordinates()
ls = 300
ws = 300
hs = 14
mesh = BoxMesh(Point(0,0,0), Point(0.014,0.012,0.000242),ls-1,ws-1,hs-1)
coor1 = mesh.coordinates()
lr = coor1[:,0].max()
ll = coor1[:,0].min()
wb = coor1[:,1].min()
wt = coor1[:,1].max()
hmax = coor1[:,2].max()
hmin = coor1[:,2].min()
# read the mesh of block7 to save coordinate
####to generate the simulation mesh
T = 1200*3.125e-6                                                          # T is final time or total time (s)
V = FunctionSpace(mesh, 'P', 1)


num_steps = 1200                                                 # num_steps is the number of time steps
t = 0
dt = T/num_steps                                                # time step size 
Rb = 2.26                                                     # the thermal resistor of convectional face
#ch = 1/(Rb*Ach)                                                         # ch is convective coefficient
h_c = 2.40598e4
Ta = 0                                                         # Ta is initial /ambient temperature 
#Nu = 5                                                          # Nu is the number if functional unit
thick_actl =  0.0000557976                                        # lthick_actl is the total thickness of active layer (m)
#import power trace file
#import floorplan file
#Line Format: <unit-name>\t<width>\t<height>\t<left-x>\t<bottom-y>\t[<specific-heat>]\t[<resistivity>]
#compute power desity 
#create mesh and define function space 

#import floorplan file
#Line Format: <unit-name>\t<width>\t<height>\t<left-x>\t<bottom-y>\t[<specific-heat>]\t[<resistivity>]
#compute power desity


#define thermal conductivity 
tol = 1E-14
k_0 = 100														# silicon conductivity      (w/(m.k))
k_1 = 100                                                   #oxide silicon thermal conductivity
kappa = Expression('x[2] <= 0.00045 + tol ? k_0 : k_1', degree=0,tol=tol, k_0=k_0, k_1=k_1) #define subdomain
#define density 
D0 = 2.33e3                                                         # silicon density    (kg/m^3)
D1 = 2.33e3                                                        # oxide silicon density
DS1 = Expression('x[2] <= 0.00045 + tol ? D0 : D1', degree=0,tol=tol, D0=D0, D1=D1) #define subdomain
#define specific heat
c0 = 751.1                                                         # silicon specific heat   (J/(kg.k))
c1 = 751.1                                                         # oxide silicon specific heat
sc = Expression('x[2] <= 0.00045 + tol ? c0 : c1', degree=0,tol=tol, c0=c0, c1=c1) #define subdomain

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
solution = []
for n in range(0,num_steps):
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
v2d = vertex_to_dof_map(V)
if N_STATE == 0:
    a = DS1*sc*u*v*dx + dt*kappa*dot(grad(u), grad(v))*dx + dt*sum(integrals_R_a)
	# Compute solution
    u = Function(V)
    dd = 0
	# define the data file name
    header_data = 'Block'+str(N_index)+'_sol_pred'
    data_file_name = header_data + '.txt'
    for n in range(0,num_steps):# changed to save time
		#update current power source term
        print(dd)
        starttime = time.time()
        solution_load_file_name ="solution_AMD_CPU/file_" + str(n) + "h5"
        solution_file = HDF5File(mesh.mpi_comm(), solution_load_file_name, "r")
        solution_file.read(u, "solution")
        data_file = open(data_file_name,'a')
        for kk in coor2:
            data_file.write('%.16g\t' % (u(kk)))
        data_file.write('\n')
        data_file.close() 
        dd += 1
        solution_file.close()
