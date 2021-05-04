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
mode = 1; # First mode number (Note: modes i and j = 0, excitation only on the "z" axis)
k_closed(n,L) = (n/2)*2*pi/L; # Resonance wavenumber, k = ω/c, where c is the speed of propagation and ω is the angular frequency. k_closed
k_cup(n,L) = (2*n-1)*pi/(2*L); #k_cup
phi_cup(k,x) = sin.(k.*x) #closed-open
phi_closed(k,x) = cos.(k.*x) #closed-closed BCs
phi_open(k,x) = sin.(k.*x) 

k = k_closed(mode,L)

CW = 343; # Speed of sound in m/s
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
BCFace = [1. 1. 0.
          2. 1. 0.
          3. 1. 1.
          4. 1. 0.
          5. 1. 0.
          6. 1. 0.]; #closed
#mshd = "./data/"
files = ["cube_coarsest.msh" "cube_coarsest.msh" "cube_coarse.msh" "cube_fine.msh" "cube_finest2.msh" "cube_finest3.msh"]
#files = ["cube_coarsest.msh" "cube_finest3.msh"]

t = [];
th=[];
Gnelem = [];
#Res = [];
#Resh = [];

for i in files[1:6]
    mesh = meshio.read(i)
    NOS_GEO = [1:size(mesh.points,1) mesh.points]
    nelem = size(mesh.cells_dict["triangle"],1)
    ELEM = [1:nelem mesh.cells_dict["triangle"].+1 mesh.cell_data["gmsh:geometrical"][1]]
    CDC,NOS = const3D_tri.gera_vars(ELEM,BCFace,NOS_GEO)
    println("Malha ",i,", nelem = ",nelem)
    elemint = []
    info = [NOS_GEO,ELEM,elemint,CDC] 
    
    #Hmatrix
#    tsolve = @elapsed TH,qH,T_pintH,qzH = const3D_tri.solveH(info,PONTOS_int,BCFace,k,true)
#    println("tsolve = ",tsolve)
    
    #Conventional
    tsolvec = @elapsed Tc,qc,T_pint,qz = const3D_tri.solve(info,PONTOS_int,BCFace,k,true)
    println("tsolvec = ",tsolvec)
    
#    global Resh = append!(Resh, [TH, T_pintH])
#    global Res = append!(Res, [Tc, T_pint])
    
    global t = append!(t,tsolvec)
#    global th = append!(th,tsolve)
    global Gnelem = append!(Gnelem,nelem)
end # files for

#save(string("./data/convergence_julia.jld"),"t",t,"th",th,"Gnelem",Gnelem)
save("convergence_julia.jld","t",t,"Gnelem",Gnelem)
