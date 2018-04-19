# Boundary element method implementation for the Helmholtz equation using constant bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Modules necessary: SpecialFunctions.jl
using SpecialFunctions
using PyPlot
include("dep.jl") # Includes the dependencies
FR = 89.7972 # Frequency of the problem [rad/s]
CW = 343 # Wave propagation speed [m/s]
k = FR/CW # Wave number
# Gaussian quadrature - generation of points and weights [-1,1]
npg=6; # Number of integration points
qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
println("Importing mesh...")
#NOS_GEO,ELEM,elemint,CDC = lermsh("../../dados/cilindro_10x100mm.msh",3) #Read the mesh generated using Gmsh
include("dad1.jl")
NOS = mostra_geo(NOS_GEO,ELEM) #Generate the physical nodes for constant elements
nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
println("Building G and H matrices...")
@time G,H,phi_inc = cal_GeH(NOS,NOS_GEO,ELEM,k,qsi,w,0) #Compute the G and H matrices
println("Applying boundary conditions to build A and b for the linear system...")
@time A,b = aplica_cdc(G,H,CDC) #Applies the boundary conditions and returns matrix A and vector b for the linear system
println("Solving the linear system...")
@time x = A\b # Solves the linear system
println("Separating acoustic pressure from flux...")
@time phi,q = monta_Teq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
T_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
