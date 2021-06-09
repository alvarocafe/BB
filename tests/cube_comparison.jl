# Boundary element method implementation for the Helmholtz and Laplace
#equations using constant  bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Contains the dependencies for the linear constant element integration.
#The main function is const2D.solve() which builds the influence matrices,
#applies the boundary conditions, solves the linear system and returns the
#value of the potential and its gradient at boundary and domain points.
using PyCall, Plots, JLD
include("../BEM_base.jl")

### Analytical solution
L = 10; # Square length
T_square(x,L) = 1 .- x./L;
k = 0
inc = [0]; # There'll be no incident wave

# BEM modelling
# Gaussian quadrature - generation of points and weights [-1,1]
npg=6; # Number of integration points
qsi,w = const3D_tri.Gauss_Legendre(-1,1,npg) # Generation of the points and weights
# Python - mesh.io
meshio = pyimport("meshio")
# Build the domain points
n_pint = 40;
PONTOS_int = zeros(n_pint,4);
delta = 0.1; # distance from both ends 
passo = (L-2*delta)/(n_pint-1);
for i = 1:n_pint
    PONTOS_int[i,:] = [i 0 0 delta+(i-1)*passo];
end
# Set the boundary conditions for each face. A cube has 6 faces
BCFace = [1. 0. 0.
          2. 1. 0.
          3. 0. 1.
          4. 1. 0.
          5. 1. 0.
          6. 1. 0.]; #closed
mshd = "./data/"
files = ["cube_coarsest.msh" "cube_coarsest.msh" "cube_coarse.msh" "cube_fine.msh" "cube_finest2.msh" "cube_finest3.msh"]

i = string(mshd,files[1])
print(i)
mesh = meshio.read(i)
NOS_GEO = [1:size(mesh.points,1) mesh.points]
nelem = size(mesh.cells_dict["triangle"],1)
ELEM = [1:nelem mesh.cells_dict["triangle"].+1 mesh.cell_data["gmsh:geometrical"][1]]
CDC,NOS = const3D_tri.gera_vars(ELEM,BCFace,NOS_GEO)
println("Malha ",i,", nelem = ",nelem)
elemint = []
info = [NOS_GEO,ELEM,elemint,CDC] 

# Helmholtz
#Hmatrix
tsolve = @elapsed pH,qpH,p_pintH,qpzH = const3D_tri.solveH(info,PONTOS_int,BCFace,k,true)
println("tsolve = ",tsolve)
#Conventional
tsolvec = @elapsed pc,qpc,p_pint,qpz = const3D_tri.solve(info,PONTOS_int,BCFace,k,true)
println("tsolvec = ",tsolvec)

# # Potencial
# #Hmatrix
# tsolve = @elapsed TH,qH,T_pintH,qzH = potconst3D.solveH(info,PONTOS_int,BCFace,k)
# println("tsolve = ",tsolve)
# #Conventional
# tsolvec = @elapsed Tc,qc,T_pint,qz = potconst3D.solve(info,PONTOS_int,BCFace,k)
# println("tsolvec = ",tsolvec)

plot(PONTOS_int[:,4],real(p_pint),label="BEM convencional")
plot!(PONTOS_int[:,4],T_square(PONTOS_int[:,4],L),label="Solução analítica")
plot!(PONTOS_int[:,4],real(p_pintH),label="Matrizes Hierárquicas")
