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
mesh_file_name = "../gmsh_mesh/Block" + str(N_index) + ".xdmf"
mesh_file = XDMFFile(comm, mesh_file_name)
mesh = Mesh()
mesh_file.read(mesh)
V = FunctionSpace(mesh, 'P', 1)
coor = mesh.coordinates()
v2d = vertex_to_dof_map(V)
lr = coor[:,0].max()
ll = coor[:,0].min()
wb = coor[:,1].min()
wt = coor[:,1].max()
hmax = coor[:,2].max()
hmin = coor[:,2].min()
h_c = 2.40598e4
Ta = 0
N_mode = 30
Nu = 13 
# Nu is the number if functional unit
N_index1 = int(sys.argv[1])   # N_index1 is the block whose gradient is computed
N_index2 = int(sys.argv[2]) # the pod mode is used

config_data = loadtxt("../Glue_block/config_block.txt",dtype="int")
if config_data[N_index1-1,N_index2+1] != 30:
    f_index1 = int(config_data[N_index1-1,N_index2+1])
    if f_index1 > 1:
        f_index2 = int(5 - f_index1)
    else:
        f_index2 = int(1 - f_index1)
else:
    exit()
flp = loadtxt('../Floorplan_AMD_multiblock_cutting.txt')

#DX = flp[4][0] - flp[8][0]
#DY = flp[8][3] - flp[4][3]
file_name1 = "../gmsh_mesh/Block_" + str(N_index1)+"_intf_"+str(f_index1)+".txt"
interface1 = loadtxt(file_name1)
file_name2 = "../gmsh_mesh/Block_" + str(N_index2)+"_intf_"+str(f_index2)+".txt"
interface2 = loadtxt(file_name2)

##################################adjust the index of N_index2 to map on the interface of N_index1
tol = 1E-14
adj_index1 = []
adj_index2 = []
for i in range (0,len(interface1)):
    for j in range (0,len(interface2)):
        if abs(interface2[j][1] - interface1[i][1] ) < tol and abs(interface2[j][2] - interface1[i][2]) < tol and abs(interface2[j][3] - interface1[i][3]) < tol:
            adj_index1.append(int(interface1[i][0]))
            adj_index2.append(int(interface2[j][0]))

#Num_nodes = mesh.num_vertices()
# ########################################################################################################################################
#                              import geometry and mesh from Gmsh                                                                        #
# ########################################################################################################################################
# import mesh file ---------multiblock.msh
#define thermal conductivity 
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
for n in range(0,N_mode):
	podmode.append('u1')
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
filename3 = "./buidling_blk_" +str(N_index2)+"/podmode_fenics.txt"
podmode2 =loadtxt(filename3)  
#loadtxt('../../Block9/pod_data/podmode_fenics_flp9.txt') # read the podmode from text file saved from HDF5 file and each row stand for each mode
#for i in range(0, len(adj_index2)):
#    u_n.vector()[v2d[adj_index2[i]]]=100
# if we have the solutiuon, we just need to  Load solution# #######################read the solution from a solution file ######################

u = Function(V)
for n in range(0,N_mode):
	podmode_load_file_name = "./buidling_blk_" + str(N_index1)+"/podmode/mode_" + str(n) + "h5"
	podmode_file = HDF5File(mesh.mpi_comm(), podmode_load_file_name, "r")
	podmode_file.read(u, "solution")
	podmode[n] = interpolate(u0,V)
	podmode[n].assign(u)
podmode_file.close()
"""
for n in range(0,N_mode):
    podmode_load_file_name = "./finite_difference_method/podmode_grad_Y/mode_" + str(n) + "h5"
    podmode_file = HDF5File(mesh.mpi_comm(), podmode_load_file_name, "r")
    podmode_file.read(u, "solution")
    grad_y[n] = interpolate(u0,V)
    grad_y[n].assign(u)
    podmode_file.close()
"""

DG = FunctionSpace(mesh, "DG", 0)
###########################compute the diagonal-Boundary term of G########################
n = FacetNormal(mesh)
#hh = CellDiameter(mesh)
G_BC_matrix = np.zeros((N_mode, N_mode))
G_BC_matrix_penalty = np.zeros((N_mode, N_mode))
for i in range (0,N_mode):
    print(i)
    u_n = interpolate(u0,V)
    for u_in in range(0, len(adj_index2)):
        u_n.vector()[v2d[adj_index1[u_in]]]=podmode2[i][adj_index2[u_in]]
    u_n_cell =interpolate(u_n,DG)
    #F = dot(u_n,inner(grad(podmode[j]),n))*ds(f_index1)
    for j in range (0,N_mode):
        #G_BC_matrix_penalty[i][j] = assemble(u_n*podmode[j]*ds(3))
        G_BC_matrix_penalty[i][j] = assemble(dot(u_n,podmode[j])*ds(f_index1))
        #print(G_BC_matrix_penalty[i][j])
        if f_index1 > 1:
            F = dot(u_n_cell,grad(podmode[j])[1])*dx
        else:
            F = dot(u_n_cell,grad(podmode[j])[0])*dx
        for cell in cells(mesh):
            M = cell.midpoint()
            if f_index1 == 0:
                if M.x()< ll + 1.19e-4:
                    G_BC_matrix[i][j] += (-assemble_local(F,cell))
            elif f_index1 == 1:
                if M.x() > lr - 1.19e-4:
                    G_BC_matrix[i][j] += (assemble_local(F,cell))
            elif f_index1 == 2:
                if M.y()< wb + 1.19e-4:
                    G_BC_matrix[i][j] += (-assemble_local(F,cell))
            else:
                if M.y() > wt - 1.19e-4:
                    G_BC_matrix[i][j] += (assemble_local(F,cell))
        #G_BC_matrix_penalty[i][j] = assemble(dot(u_n,podmode[j])*ds(f_index1))
        #G_BC_matrix_penalty[i][j] = assemble(dot(u_n,podmode[j]/hh)*ds(f_index1))

##############save the G_matrix_BC_diag
header = './buidling_blk_' + str(N_index1)+'/G_matrix_offdiag' + '_' + str(N_index1) +'_'+str(N_index2)
#header = 'G_matrix_offdiag' + str(N_index2)+'_' + str(N_index1)
G_BC_matrix_file_name = header + '.txt'
G_BC_matrix_file = open(G_BC_matrix_file_name,'w')
for i in range(0,N_mode):
    for j in range(0, N_mode):
        G_BC_matrix_file.write('%.16g\t' % (G_BC_matrix[i][j]))
    G_BC_matrix_file.write('%\n')
G_BC_matrix_file.close()

##############save the G_matrix_BC_diag
header = './buidling_blk_' + str(N_index1)+'/G_matrix_offdiag_penalty'+'_' + str(N_index1) +'_'+str(N_index2)
#header = 'G_matrix_offdiag_penalty'+ str(N_index2)+'_' + str(N_index1)
G_BC_matrix_penalty_file_name = header + '.txt'
G_BC_matrix_penalty_file = open(G_BC_matrix_penalty_file_name,'w')
for i in range(0,N_mode):
    for j in range(0, N_mode):
        G_BC_matrix_penalty_file.write('%.16g\t' % (G_BC_matrix_penalty[i][j]))
    G_BC_matrix_penalty_file.write('%\n')
G_BC_matrix_penalty_file.close()

#vtkfile = File("data_assign.pvd")
#vtkfile << u_n



"""
import matplotlib.pyplot as plt
P = plot(u_n)
plt.colorbar(P)
#plot(mesh)
plt.show()
"""
