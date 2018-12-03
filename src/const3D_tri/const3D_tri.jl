# Boundary element method implementation for the Helmholtz equation using constant tridimensional triangular elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Module for the constant three-dimensional triangular element
# Contains the dependencies for the triangular element integration. The main function is const3D_tri.solve() which builds the influence matrices, applies the boundary conditions, solves the linear system and returns the value of the velocity potential and its flux at boundary and domain points.

module const3D_tri
using PyPlot
using PyCall
# @pyimport matplotlib.colors as col
@pyimport matplotlib.cm as cm
@pyimport mpl_toolkits.mplot3d as mp
@pyimport mpl_toolkits.mplot3d.art3d as ar
plt=PyPlot

include("format.jl") # curve interpolation formatting
include("cal.jl") # element integration calculating functions
include("H_mat.jl") # H-Matrices support for building the cluster tree and blocks
include("interp.jl") # approximation  using Lagrange polynomial interpolation
include("ACA.jl") # approximation using ACA


function solve(info,PONTOS_int,BCFace,k)
    NOS_GEO,ELEM,elemint,CDC = info
    NOS = mostra_geoTRI(NOS_GEO,ELEM) #Generate the physical nodes for constant elements
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    # If there's a defined BCFace, the boundary condition matrix is created. If no BCFace is defined, then the mesh reader has already built the boudnary conditions matrix.
    if isempty(BCFace) != true
	CDC = gera_CDC(ELEM,BCFace); #Monta a matriz de condicoes de contorno
    end
    # Gaussian quadrature - generation of points and weights [-1,1]
    npg=12; # Number of integration points
    qsi,w = Gauss_Legendre(0,1,npg) # Generation of the points and weights
    println("Building G and H matrices...")
    @time G,H,phi_inc = cal_GeH(NOS,NOS_GEO,ELEM,k,qsi,w,0) #Compute the G and H matrices
    println("Applying boundary conditions to build A and b for the linear system...")
    @time A,b = aplica_cdc(G,H,CDC) #Applies the boundary conditions and returns matrix A and vector b for the linear system
    println("Solving the linear system...")
    @time x = A\b # Solves the linear system
    println("Separating acoustic pressure from flux...")
    @time phi,q = monta_Teq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
    println("Solving for domain points.")
    @time phi_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
    @time dphidx,dphidy,dphidz=calc_q_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
    return phi,q,phi_pint,dphidz
end

function H_solve(info,PONTOS_int,BCFace,k)
    NOS_GEO,ELEM,elemint,CDC = info
    NOS = mostra_geoTRI(NOS_GEO,ELEM) #Generate the physical nodes for constant elements
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    # If there's a defined BCFace, the boundary condition matrix is created. If no BCFace is defined, then the mesh reader has already built the boudnary conditions matrix.
    if isempty(BCFace) != true
	CDC = gera_CDC(ELEM,BCFace); #Monta a matriz de condicoes de contorno
    end
    # Gaussian quadrature - generation of points and weights [-1,1]
    npg=12; # Number of integration points
    qsi,w = Gauss_Legendre(0,1,npg) # Generation of the points and weights
    A,b = interp()
    println("Separating acoustic pressure from flux...")
    @time phi,q = monta_Teq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
    println("Solving for domain points.")
    @time phi_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
    @time dphidx,dphidy,dphidz=calc_q_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
    return phi,q,phi_pint,dphidz
end
end
