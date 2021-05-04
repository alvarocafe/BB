# Boundary element method implementation for the Helmholtz and Laplace
#equations using constant  bidimensional elements
# Authors: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com,
# Fernando Barreto Soares - fernando.bsoares3@gmail.com
# Contains the dependencies for the linear constant element integration.
#The main function is const2D.solve() which builds the influence matrices,
#applies the boundary conditions, solves the linear system and returns the
#value of the potential and its gradient at boundary and domain points.
using PyCall, Plots, JLD
include("../BEM_base.jl")

### Analytical solution
L = 10; # Square length
mode = 2; # First mode number (Note: modes i and j = 0, excitation only on the "z" axis)
k_closed(n,L) = (n/2)*2*pi/L; # Resonance wavenumber, k = ω/c, where c is the speed of propagation and ω is the angular frequency. k_closed
k_cup(n,L) = (2*n-1)*pi/(2*L); #k_cup
phi_cup(k,x) = sin.(k.*x) #closed-open
phi_closed(k,x) = cos.(k.*x) #closed-closed BCs
phi_open(k,x) = sin.(k.*x)

CW = 343; # Speed of sound in m/s
inc = [0]; # There'll be no incident wave

# BEM modelling
# Gaussian quadrature - generation of points and weights [-1,1]
npg=6; # Number of integration points
qsi,w = const3D_tri.Gauss_Legendre(-1,1,npg) # Generation of the points and weights
# Python - mesh.io
meshio = pyimport("meshio")
# Build the domain points
L = 10; # Length of the cube
n_pint = 40;
PONTOS_int = zeros(n_pint,4);
delta = 0.1; # distance from both ends 
passo = (L-2*delta)/(n_pint-1);
for i = 1:n_pint
    PONTOS_int[i,:] = [i 0 0 delta+(i-1)*passo];
end
# Set the boundary conditions for each face. Faces 1 and 3 are perpendicular do "z" axis, 4 and 2 are perpendicular to "x" and 6 and 3 are perpendicular to "y".
#BCFace = [1. 1. 0.
#          2. 1. 0.
#          3. 0. 1.
#          4. 1. 0.
#          5. 1. 0.
#          6. 1. 0.]; #cup
#BCFace = [1. 0. -1.
#          2. 1. 0.
#          3. 0. 1.
#          4. 1. 0.
#          5. 1. 0.
#          6. 1. 0.];
BCFace = [1. 0. 0.
          2. 1. 0.
          3. 0. 1.
          4. 1. 0.
          5. 1. 0.
          6. 1. 0.]; # open
#BCFace = [1. 1. -1.
#          2. 1. 0.
#          3. 1. 1.
#          4. 1. 0.
#          5. 1. 0.
#          6. 1. 0.]; 
#BCFace = [1. 1. 0.
#          2. 1. 0.
#          3. 1. 1.
#          4. 1. 0.
#          5. 1. 0.
#          6. 1. 0.];  #closed_closed
mshd = "./data/"
file = "cube_coarse.msh"
#t = [];
#Gnelem = [];
mesh = meshio.read(string(mshd,file))
NOS_GEO = [1:size(mesh.points,1) mesh.points]
nelem = size(mesh.cells_dict["triangle"],1)
ELEM = [1:nelem mesh.cells_dict["triangle"].+1 mesh.cell_data["gmsh:geometrical"][1]]
CDC,NOS = const3D_tri.gera_vars(ELEM,BCFace,NOS_GEO)
println("Malha ",file,", nelem = ",nelem)
elemint = []
info = [NOS_GEO,ELEM,elemint,CDC]

#ks = [k_closed(1,L), k_closed(2,L), k_closed(3,L)]
#T_pintc = zeros(ComplexF64,n_pint,length(ks))
#T_analitic = zeros(n_pint,length(ks))

#for i = 1:length(ks)
#  k = ks[i]
#  Tc,qc,T_pint,qz = const3D_tri.solve(info,PONTOS_int,BCFace,k,true)
#  T_pintc[:,i] = T_pint
#  T_analitic[:,i] = phi_open(k,PONTOS_int[:,4])
#end

k = k_closed(mode,L)

Tc,qc,T_pint,qz = const3D_tri.solve(info,PONTOS_int,BCFace,k,true)

T_analitic = phi_open(k_closed(mode,L),PONTOS_int[:,4])

#errorms = zeros(length(ks))
#for i = 1:length(ks)
#  errorms[i] = sqrt(sum((real(T_pintc[:,i])/maximum(real(T_pintc[:,i]))-T_analitic[:,i]).^2)/n_pint)/(maximum(T_analitic[:,i])-minimum(T_analitic[:,i])) #pg.74 dissertação lucas campos
#end  	

#save(string("./data/assymetric_open.jld"),"T_pintc",T_pintc)
#save(string("./data/errorms_open.jld"),"errorms",errorms)

#TH,qH,T_pintH,qzH = const3D_tri.solveH(info,PONTOS_int,BCFace,k,true)

#Results plotting

pyplot()

#use the regular plot to obtain the pressure distribution inside the cube at a given frequency (Line Plot)
plot(PONTOS_int[:,4],real(T_pint)/maximum(real(T_pint)),xaxis=("z"),yaxis=("u"),title=(string("Solução numérica de u, n = ",mode)),label=("MEC"))
plot!(PONTOS_int[:,4],T_analitic,label=("Analítico"))
#plot(PONTOS_int[:,4],real(T_pintH)/maximum(real(T_pintH)),xaxis=("z"),yaxis=("u"),title=(string("Solução numérica de u, n = ",mode)),label=("MEC Hmatrix"))
