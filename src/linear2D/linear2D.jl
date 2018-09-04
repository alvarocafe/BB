# Boundary element method implementation for the Helmholtz and Laplace equations using linear  bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Contains the dependencies for the linear element integration.
#The main function is linear2D.solve() which builds the influence matrices, applies the boundary conditions,
#solves the linear system and returns the value of the velocity potential and its flux at boundary and domain points.

module linear2D
using SpecialFunctions
using PyPlot

include("dep.jl") # Includes the dependencies


function solve()
    G,H,q = monta_GeH(ELEM_d,NOS_d,CDC,k,fc,qsi,w); # Evaluates the influence matrices H and G and source influence array q
    A,b=aplica_CDC(G,H,CDC);	# Applies the boundary conditions to build the linear system matrix A and vector b
    x=A\(b); # Solves the linear system
    phi,qphi = Monta_Teq(x,CDC);	# Applies the boundary conditions to build the velocity potential and its flux
    return phi, qphi    
end

end # end module linear2D
