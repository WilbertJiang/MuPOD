# Welcome to MuPOD!
MuPOD is a data-driven thermal simulator that can be dynamically reorganized to adapt to variations in chip floorplans.
This `README.md` provides step-by-step guidance on installing and using MuPOD. The simulator can also be seamlessly integrated into EDA toolchains as an efficient and accurate thermal model.
This work was recognized with the Prof. Avram Bar-Cohen [Best Paper Award at ITherm 2022](https://www.ieee-itherm.net/2022-best-paper-winners/).

<p align="center">
  <img src="/Image/BPA.jpg" alt="Best Paper Award" width="600">
</pr>

If you use any component of MuPOD in your work, please cite:

```
[1] L. Jiang, A. Dowling, Y. Liu and M. -C. Cheng, "Chip-level thermal simulation for a multicore processor using a multi-block model enabled by proper orthogonal decomposition"
ITherm (2022), p. 2022
```


# Overview
MuPOD combines Proper Orthogonal Decomposition (POD) with domain decomposition:

**1. Partitioning into Blocks**:
The chip is divided into smaller building blocks (cores, caches, I/O, memory, etc.) based on its floorplan. Each block is simulated individually with a high-resolution FEM tool [FEniCS](https://fenicsproject.org/).

**2. Training POD Modes**:
For each block, dynamic thermal data is collected under varying power and boundary conditions. POD extracts a small set of basis functions (modes) capturing most of the thermal behavior. These modes form a reduced-order model (ROM) for the block.

**3. Dynamically Assembling the Multi-Block Model**:
Block-level POD models are combined into a chip-level simulator for the entire chip. At block interfaces, the discontinuous Galerkin (DG) method enforces thermal continuity by balancing temperature and heat flux across boundaries.

To download and install MuPOD
```
git clone --recursive https://github.com/WilbertJiang/MuPOD.git
```
# How to Install MuPOD and Run MuPOD Step by Step
## 1. Dependencies
MuPOD is developed on the FEniCS platform, which provides a flexible framework for solving partial differential equations (PDEs) using finite element methods. FEniCS should be pre-installed using the following command:  
### FEniCS install
```
sudo apt-get install --no-install-recommends software-properties-common  
sudo add-apt-repository ppa:fenics-packages/fenics  
sudo apt-get update  
sudo apt-get install fenics
```
Please refer to the FEniCS installation guide for more detailed instructions on installation and troubleshooting: [FEniCS download](https://fenicsproject.org/download/.).
### Building tools installation   
To run the C++ version FEniCS, you need to make sure that the build tools are installed
```
sudo apt install cmake make g++ -y
```
### C++ FEniCS installation
If the cmake are installed on your server, you can then run the following commands to install C++ version FEniCS
```
sudo apt-get install --no-install-recommends software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install --no-install-recommends fenics
```

## 2. Training data collection
The temperature data required for training an individual POD model in MuPOD can be obtained through either experimental measurements or numerical simulations. For example, the finite element method (FEM) implemented in  [FEniCS](https://fenicsproject.org/), which is used in [PODTherm-GP](https://github.com/WilbertJiang/PODTherm_GP), can also be employed to generate the training data for MuPOD. 

This can be done by navigating to the home folder of PODTherm-GP and running the following command:
```
    cd ./src  
    ffc -l dolfin Space.ufl  
    cd ..  
    mkdir build  
    cd ./build  
    cmake ..  
    make 
 ```
 Then the executable file can be run with one or multiple processes. For instance, the component of training data collection can be performed by 
 ```
 mpirun -n 20 ./Therm_FEM
 ```
where 20 is the number of processes. 

Or, you can just run the following command to collect data
```
mpirun -n 10 python3 Sol_simu_block.py
```
## 3. Individual POD Model Training
To train an individual POD model for each building block, run 
```
mpirun -n 10 python3 ComputingAM_flp7.py
python3 save_podmode.py
```
Or, just run 
```
./generate_mode.sh
```
Once this command is done, an individual POD model is generated for each building block, including the thermal conductance matrix: **G**, thermal capacitance matrix: **C**, and power density in POD space: **P**
## 4. Dynamically Assembling
With the generated individual POD models, a thermal model can be dynamically constructed for the entire chip
```
./run.sh
```
Using this command, 
And then, the POD models can be dynamically assembled with the variations of floorplan 
"./Construct/run2.sh"


