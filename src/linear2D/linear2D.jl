# Boundary element method implementation for the Helmholtz and Laplace
#equations using linear discontinuous bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Contains the dependencies for the linear discontinuous element integration.
#The main function is desclin2D.solve() which builds the influence matrices,
#applies the boundary conditions, solves the linear system and returns the
#value of the potential and its gradient at boundary and domain points.

module desclin2D
using SpecialFunctions
using KrylovMethods
using PyCall
using PyPlot
plt=PyPlot
@pyimport matplotlib.tri as tri


include("format.jl") # curve interpolation formatting
include("kernel.jl") # 
include("cal.jl") #
include("H_mat.jl") # H-Matrices support for building the cluster tree and blocks
include("interp.jl") # H-Matrices using Lagrange polynomial interpolation
include("ACA.jl") # H-Matrices using ACA


function solve(info,PONTOS_int,fc,BCFace,k)
## CBIE - Conventional Boundary Integral Equation
	NOS_GEO,NOS,ELEM,CDC = info;
	nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
	b1 = 1:nnos # Array containing all the indexes for nodes and elements which will be used for integration
	# Gaussian quadrature - generation of points and weights [-1,1]
	npg=6; # Number of integration points
	qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
#	println("Building A and b matrices using the traditional colocation BEM for constant elements.")
	@time A,b = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])  # Builds A and B matrices using the collocation technique and applying the boundary conditions
	x = A\b # Solves the linear system
	phi,qphi = monta_phieq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
#	println("Evaluating at domain points.")
	@time phi_pint = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,fc,qsi,w,k) # Evaluates the value at domain points
return phi, qphi, phi_pint, phi_pint
end

function solveH(info,PONTOS_int,fc,BCFace,k)
## H-Matrix BEM - Interpolation using Lagrange polynomial
	NOS,NOS_GEO,ELEM,CDC = info;
#	println("Building Tree and blocks using H-Matrices.")
	@time Tree,block = cluster(NOS[:,2:3],floor(sqrt(length(NOS))),2)
	# Gaussian quadrature - generation of points and weights [-1,1]
	npg=6; # Number of integration points
	qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
#	println("Building A and b matrices using H-Matrix with interpolation.")
	@time Ai,bi = Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
	xi = gmres(vet->matvec(Ai,vet,block,Tree),bi,5,tol=1e-5,maxIter=1000,out=0) #GMRES nas matrizes do ACA
	phii,qphii = monta_phieq(CDC,xi[1]) # Applies the boundary conditions to return the velocity potential and flux
#	println("Evaluating values at internal points.")
	@time phi_pinti = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phii,qphii,fc,fc,qsi,w,k) # Evaluates the value at internal (or external) points
return phii,qphii,phi_pinti,phi_pinti
end

function solvepot(info,PONTOS_int,fc,BCFace,k)
## CBIE - Conventional Boundary Integral Equation
	NOS_GEO,NOS,ELEM,CDC = info;
	nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
	b1 = 1:nnos # Array containing all the indexes for nodes and elements which will be used for integration
	# Gaussian quadrature - generation of points and weights [-1,1]
	npg=6; # Number of integration points
	qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
#	println("Building A and b matrices using the traditional colocation BEM for constant elements.")
	@time A,b = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])  # Builds A and B matrices using the collocation technique and applying the boundary conditions
	x = A\b # Solves the linear system
	phi,qphi = monta_phieq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
#	println("Evaluating at domain points.")
	@time phi_pint = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,fc,qsi,w,k) # Evaluates the value at domain points
return phi, qphi, phi_pint, phi_pint
end

function solveHpot(info,PONTOS_int,fc,BCFace,k)
## H-Matrix BEM - Interpolation using Lagrange polynomial
	NOS,NOS_GEO,ELEM,CDC = info;
#	println("Building Tree and blocks using H-Matrices.")
	@time Tree,block = cluster(NOS[:,2:3],floor(sqrt(length(NOS))),2)
	# Gaussian quadrature - generation of points and weights [-1,1]
	npg=6; # Number of integration points
	qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
#	println("Building A and b matrices using H-Matrix with interpolation.")
	@time Ai,bi = Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
	xi = gmres(vet->matvec(Ai,vet,block,Tree),bi,5,tol=1e-5,maxIter=1000,out=0) #GMRES nas matrizes do ACA
	phii,qphii = monta_phieq(CDC,xi[1]) # Applies the boundary conditions to return the velocity potential and flux
#	println("Evaluating values at internal points.")
	@time phi_pinti = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phii,qphii,fc,fc,qsi,w,k) # Evaluates the value at internal (or external) points
return phii,qphii,phi_pinti,phi_pinti
end

end
