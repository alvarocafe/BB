# Boundary element method implementation for the Helmholtz and Laplace equations using linear discontinuous  bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Contains the dependencies for the linear discontinuous element integration.
#The main function is desclin2D.solve() which builds the influence matrices, applies the boundary conditions,
#solves the linear system and returns the value of the velocity potential and its flux at boundary and domain points.

module desclin2D
using SpecialFunctions
using PyPlot

include("dep.jl") # Includes the dependencies
include("H_mat.jl") # H-Matrices support for building the cluster tree and blocks
include("beminterp.jl") # H-Matrices using Lagrange polynomial interpolation
include("ACA.jl") # H-Matrices using ACA

function solve(info,PONTOS_int,fc,BCFace,k)
    ## CBIE - Conventional Boundary Integral Equation
    NOS_GEO,NOS,ELEM,CDC = info;
    ne = size(ELEM,1);
    # Gaussian quadrature - generation of points and weights [-1,1]
    npg=6; # Number of integration points
    qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    println("Building A and b matrices using the traditional colocation BEM for constant elements.")    
    # Evaluates the continuous shape functions at the points where the descontinuous nodes will be created. This nodes will be used with the descontinuous linear shape function later on. 
    N11,N21=calc_fforma(-1/2); 
    N12,N22=calc_fforma(1/2);
    # Allocates the new nodes, connectivity and boundary conditions matrices 
    NOS_d = zeros(2*ne,3);
    ELEM_d = zeros(Int,ne,3);
    CDC_d = zeros(2*ne,3);
    j = 0;
    for i=1:ne # Loop over the elements
        no1::Int = ELEM[i,2];	# First point of the element
        no2::Int = ELEM[i,3]	# Second point of the element
        # Coordinate of the continuous element points
        x1 = NOS[no1,2];    y1 = NOS[no1,3];
        x2 = NOS[no2,2];    y2 = NOS[no2,3];
        j=j+1;
        NOS_d[j,:] = [j N11*x1+N21*x2 N11*y1+N21*y2]; # Evaluates the position of the first new descontinuous node
        j=j+1;
        NOS_d[j,:] = [j N12*x1+N22*x2 N12*y1+N22*y2]; # Evaluates the position of the second new descontinuous node
        ELEM_d[i,:] = [i j-1 j]; # Updates the mesh connectivity matrix
    end
    println("Evaluating matrices H, G and array q");
    G,H,q = monta_GeH(ELEM_d,NOS_d,CDC,k,fc,qsi,w); # Evaluates the influence matrices H and G and source influence array q
    A,b=aplica_CDC(G,H,CDC);	# Applies the boundary conditions to build the linear system matrix A and vector b
    println("Solving the linear system");
    x=A\(b); # Solves the linear system
    println("Reordering to obtain phi and qphi");
    phi,qphix = Monta_Teq(x,CDC);	# Applies the boundary conditions to build the velocity potential and its flux
    return phi, qphix, qphix, phi_pint 
end # end function solve

end # end module desclin2D
