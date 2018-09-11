# Boundary element method implementation for the Helmholtz and Laplace
#equations using linear discontinuous bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Contains the dependencies for the linear discontinuous element integration.
#The main function is desclin2D.solve() which builds the influence matrices,
#applies the boundary conditions, solves the linear system and returns the
#value of the potential and its gradient at boundary and domain points.
#Start----Problem--------------Method-----------------post-processing
#-----------^You are here!-------------------------------------------

module desclin2D
using SpecialFunctions, KrylovMethods, PyCall, PyPlot
plt=PyPlot
@pyimport matplotlib.tri as tri
include("format.jl") # curve interpolation formatting
include("kernel.jl") # problem statement kernel
include("cal.jl") # element integration calculating functions
include("H_mat.jl") # H-Matrices support for building the cluster tree and blocks
include("interp.jl") # approximation  using Lagrange polynomial interpolation
include("ACA.jl") # approximation using ACA

function solve(info,PONTOS_int,fc,BCFace,k)
    ## CBIE - Conventional Boundary Integral Equation
    NOS_GEO,NOS,ELEM,CDC = info
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using
    #constant elements
    b1 = 1:nnos # Array containing all the indexes for nodes and elements which
    #will be used for integration
    # Gaussian quadrature - generation of points and weights [-1,1]
    npg=6
    qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    A,b = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    x = A\b # Solves the linear system
    phi,qphi = monta_phieq(CDC,x) # Applies the boundary conditions
    phi_dom = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,fc,qsi,w,k)
    return phi, qphi, phi_domt, phi_dom
end

function solveH(info,PONTOS_int,fc,BCFace,k)
    ## H-Matrix BEM - Interpolation using Lagrange polynomial
    NOS,NOS_GEO,ELEM,CDC = info
    Tree,block = cluster(NOS[:,2:3],floor(sqrt(length(NOS))),2)
    # Gaussian quadrature - generation of points and weights [-1,1]
    npg=6
    qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    A,b = Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    x = gmres(vet->matvec(Ai,vet,block,Tree),bi,5,tol=1e-5,maxIter=1000,out=0)
    phi,qphi = monta_phieq(CDC,xi[1]) # Applies the boundary conditions
    phi_dom = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phii,qphii,fc,fc,qsi,w,k)
    return phi,qphi,phi_dom,phi_dom
end

function solvepot(info,PONTOS_int,fc,BCFace,k)
    ## CBIE - Conventional Boundary Integral Equation
    NOS_GEO,NOS,ELEM,CDC = info
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using
    #constant elements
    b1 = 1:nnos # Array containing all the indexes for nodes and elements which
    #will be used for integration
    # Gaussian quadrature - generation of points and weights [-1,1]
    npg=6; # Number of integration points
    qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    A,b = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    x = A\b # Solves the linear system
    phi,qphi = monta_phieq(CDC,x) # Applies the boundary conditions
    phi_dom = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,fc,qsi,w,k)
    return phi, qphi, phi_dom, phi_dom
end

function solveHpot(info,PONTOS_int,fc,BCFace,k)
    ## H-Matrix BEM - Interpolation using Lagrange polynomial
    NOS,NOS_GEO,ELEM,CDC = info;
    Tree,block = cluster(NOS[:,2:3],floor(sqrt(length(NOS))),2)
    # Gaussian quadrature - generation of points and weights [-1,1]
    npg=6
    qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    A,b = Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    x = gmres(vet->matvec(Ai,vet,block,Tree),bi,5,tol=1e-5,maxIter=1000,out=0)
    phi,qphi = monta_phieq(CDC,xi[1]) # Applies the boundary conditions
    phi_dom = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phii,qphii,fc,fc,qsi,w,k)
    return phi,qphi,phi_dom,phi_dom
end

end # end module desclin2D
