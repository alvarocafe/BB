# Boundary element method implementation for the Helmholtz and Laplace
#equations using constant  bidimensional elements
# Authors: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com,
# Fernando Barreto Soares - fernando.bsoares3@gmail.com
# Contains the dependencies for the linear constant element integration.
#The main function is const2D.solve() which builds the influence matrices,
#applies the boundary conditions, solves the linear system and returns the
#value of the potential and its gradient at boundary and domain points.
using JLD
using PyPlot
using PyCall
@pyimport matplotlib.colors as col
@pyimport matplotlib.cm as cm
@pyimport mpl_toolkits.mplot3d as mp
@pyimport mpl_toolkits.mplot3d.art3d as ar
#plt=PyPlot
#plot=Plots

include("../BEM_base.jl")
include("../src/const3D_tri/mostra_resultados2.jl")

freq = 40
omega = 2 * pi * freq
CW = 344
k = omega/CW

inc = [0]; # There'll be no incident wave
#The analytical solution for the pulsating sphere is
phi_sphere(k,r,a,ρ,c) = a/r*ρ*c*(-complex(0,1)*k*a/(1+complex(0,1)*k*a))*exp(complex(0,1)*k*(r-a));

# BEM modelling
# Gaussian quadrature - generation of points and weights [-1,1]
npg=6; # Number of integration points
qsi,w = const3D_tri.Gauss_Legendre(-1,1,npg) # Generation of the points and weights
# Python - mesh.io
meshio = pyimport("meshio")
# Build the domain points
L = 100; # Length of the cube
n_pint = 100;
PONTOS_int = zeros(n_pint*n_pint,4);
passo = L/(n_pint-1);
starty = -50
startz = -50
iter = 0;
for i in 1:n_pint
    for j in 1:n_pint
        global iter += 1
        PONTOS_int[iter,:] = [iter 0 starty-passo+i*passo startz-passo+j*passo]
    end
end

#BCFace = [1. 1. 0.
#          2. 1. 1.
#          3. 1. 0.
#          4. 1. -1.]; 

mshd = "./data/"
file = "sphere.msh"
#t = [];
#Gnelem = [];
mesh = meshio.read(string(mshd,file))
NOS_GEO = [1:size(mesh.points,1) mesh.points]
nelem = size(mesh.cells_dict["triangle"],1)
ELEM = [1:nelem mesh.cells_dict["triangle"].+1 mesh.cell_data["gmsh:geometrical"][1]]

geo_to_physic = [mesh.cell_data["gmsh:physical"][1] mesh.cell_data["gmsh:geometrical"][1]]

BCFace = [1 1 1] #There are 432 surfaces

CDC,NOS = const3D_tri.gera_vars(ELEM,BCFace,NOS_GEO)
println("Malha ",file,", nelem = ",nelem)
elemint = []
info = [NOS_GEO,ELEM,elemint,CDC]

#mostra_resultados2(NOS_GEO,ELEM,CDC[:,3])

Tc,qc,T_pint,qz = const3D_tri.solve(info,PONTOS_int,BCFace,k,true)

save(string("tinylev.jld"),"T_pint",T_pint)

# Results plotting

pyplot = PyPlot.plot

using ColorSchemes
pontosy = zeros(n_pint,2);
for i = 1:n_pint
    pontosy[i,:] = [i starty-passo+i*passo];
end
pontosz = zeros(n_pint,2);
for i = 1:n_pint
    pontosz[i,:] = [i startz-passo+i*passo];
end

phi_real=real(T_pint);
phi_contour = zeros(n_pint,n_pint);
phi_anal = zeros(n_pint,n_pint);
for j = 1:n_pint
    index1=n_pint*(j-1)+1
    index2=n_pint*j
    phi_contour[:,j] = phi_real[index1:index2]
    phi_anal[:,j] = phi_real[index1:index2]
end

# solar = ColorSchemes.solar.colors
# plot.contour(pontosy[:,2],pontosz[:,2],phi_contour/maximum(phi_contour),fill=true,colorbar=true,xaxis=("x"),yaxis=("z"),title=("Solução numérica TinyLev, plano y=0"),seriescolor=cgrad(ColorSchemes.viridis.colors),aspect_ratio=:equal,levels=20)

# plot.show()
#save(string("./data/sphereHnegativo.jld"),"T_pint",T_pint,"phi_contour",phi_contour)

contour(pontosy[:,2],pontosz[:,2],phi_contour/maximum(phi_contour),fill=true,colorbar=true,xaxis=("x"),yaxis=("z"),title=("Solução numérica TinyLev, plano y=0"),aspect_ratio=:equal,levels=20)

contour(pontosy[:,2],pontosz[:,2],phi_anal/maximum(phi_anal),fill=true,colorbar=true,xaxis=("x"),yaxis=("z"),title=("Solução numérica TinyLev, plano y=0"),aspect_ratio=:equal,levels=20)
