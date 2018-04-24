# Boundary element method implementation for the Helmholtz equation using constant bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Module for the constant three-dimensional quadrilateral element
# Contains the dependencies for the quadrilateral element integration. The main function is const3D_quad.solve() which builds the influence matrices, applies the boundary conditions, solves the linear system and returns the value of the velocity potential and its flux at boundary and domain points.
# Necessary Modules: SpecialFunctions.jl and KrylovMethods
using SpecialFunctions
import KrylovMethods
module const2D

include("dep.jl") # Includes the dependencies
include("beminterp.jl") # H-Matrices using Lagrange polynomial interpolation
include("ACA.jl") # H-Matrices using ACA


function solve(info,PONTOS_int,fc,BCFace,k)
## CBIE - Conventional Boundary Integral Equation
	NOS_GEO,NOS,ELEM,CDC = info;
	nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
	b1 = 1:nnos # Array containing all the indexes for nodes and elements which will be used for integration
	# Gaussian quadrature - generation of points and weights [-1,1]
	npg=6; # Number of integration points
	qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
	println("Building A and b matrices using the traditional colocation BEM for constant elements.")
	@time A,b = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])  # Builds A and B matrices using the collocation technique and applying the boundary conditions
	x = A\b # Solves the linear system
	phi,qphi = monta_phieq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
	println("Evaluating values at domain points.")
	@time phi_pint = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,fc,qsi,w,k) # Evaluates the value at domain points
return phi, qphi, phi_pint, phi_pint
end

function solveH(info,PONTOS_int,fc,BCFace,k)
import KrylovMethods
## H-Matrix BEM - Interpolation using Lagrange polynomial
	NOS,NOS_GEO,ELEM,CDC = info;
	println("Building Tree and blocks using H-Matrices.")
	@time Tree,block = cluster(NOS[:,2:3],floor(sqrt(length(NOS))),2)
	# Gaussian quadrature - generation of points and weights [-1,1]
	npg=6; # Number of integration points
	qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
	println("Building A and b matrices using H-Matrix with interpolation.")
	@time Ai,bi = Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
	xi = KrylovMethods.gmres(vet->matvec(Ai,vet,block,Tree),bi,5,tol=1e-5,maxIter=1000,out=0) #GMRES nas matrizes do ACA
	phii,qphii = monta_phieq(CDC,xi[1]) # Applies the boundary conditions to return the velocity potential and flux
	println("Evaluating values at internal points.")
	@time phi_pinti = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phii,qphii,fc,fc,qsi,w,k) # Evaluates the value at internal (or external) points
return phii,qphii,phi_pinti,phi_pinti
end
end
