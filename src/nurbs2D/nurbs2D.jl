# Boundary element method implementation for the Laplace equation using NURBS bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Necessary Modules: SpecialFunctions.jl

module nurbs2D
using SpecialFunctions
using PyPlot
include("dep.jl") # Includes the dependencies
include("H_mat.jl") # H-Matrices support for building the cluster tree and blocks
include("beminterp.jl") # H-Matrices using Lagrange polynomial interpolation
include("ACA.jl") # H-Matrices using ACA

function solve(info,PONTOS_int,fc,k)
# CBIE - conventional boundary integral equation
	collocCoord,nnos,crv,dcrv,CDC,E = info;
	#Building the problems matrices
	H, G = calcula_iso(collocCoord,nnos,crv,dcrv,E,k) # Influence matrices
	A,b= aplica_CDCiso(G,H,CDC,E);	# Applying boundary conditions
	x=A\b; # Evaluating unknown values
	Tc,qc=monta_Teqiso(CDC,x); # Separating temperature from flux
	# Applying NURBS basis functions to the values of temperature and flux
	T=E*Tc;
	q=E*qc;
	# Domain points
	Hp,Gp,phi_pint = calc_phi_pint_nurbs(PONTOS_int,collocCoord,nnos,crv,dcrv,k,Tc,qc);

end # end function solve
end # end module nurbs2D
