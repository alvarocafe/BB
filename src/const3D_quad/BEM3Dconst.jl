# Boundary element method implementation for the Helmholtz equation using constant bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Modules necessary: SpecialFunctions.jl
using SpecialFunctions
using PyPlot
include("dep.jl") # Includes the dependencies
FR = 20 # Frequency of the problem [Hz]
CW = 343*1000 # Wave propagation speed [mm/s]
k = FR/CW # Wave number
# Gaussian quadrature - generation of points and weights [-1,1]
npg=16; # Number of integration points
qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
println("Importing mesh...")
NOS_GEO,ELEM,elemint,CDC = lermsh("../../dados/cilindro_5x10mm.msh",3) #Read the mesh generated using Gmsh
NOS = mostra_geo(NOS_GEO,ELEM) #Generate the physical nodes for constant elements
nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
CCFace = [1 1 1
	  2 0 1
	  3 0 1]
CDC = gera_CDC(ELEM,CCFace); #Monta a matriz de condicoes de contorno
n_pint = 100
PONTOS_int = zeros(n_pint,4)
for i = 1:n_pint
	PONTOS_int = [i (10+(20/n_pint)*i) 0 50]
end
println("Building G and H matrices...")
@time G,H,phi_inc = cal_GeH(NOS,NOS_GEO,ELEM,FR,CW,qsi,w,0) #Compute the G and H matrices
println("Applying boundary conditions to build A and b for the linear system...")
@time A,b = aplica_cdc(G,H,CDC) #Applies the boundary conditions and returns matrix A and vector b for the linear system
println("Solving the linear system...")
@time x = A\b # Solves the linear system
println("Separating acoustic pressure from flux...")
@time phi,q = monta_Teq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
T_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,CW,FR,qsi,w,0)
