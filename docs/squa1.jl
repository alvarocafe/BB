# Boundary element method implementation for the Helmholtz equation using constant bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
include("../src/const2D/dep.jl") # Includes the dependencies
using SpecialFunctions

PONTOS = [1 0 0
	  2 6 0
	  3 6 6
	  4 0 6];
SEGMENTOS = [1 1 2 0
	     2 2 3 0
	     3 3 4 0
	     4 4 1 0];
ne = 3;
MALHA = [1 ne
	 2 ne
	 3 ne
	 4 ne];
CCSeg = [1 0 0
	 2 1 1
	 3 1 100
	 4 1 0];
PONTOS_int = [1 3 3]
NOS_GEO,NOS,ELEM,CDC = format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg) # Apply the discretization technique and builds the problems matrices for the geometrical points, physical nodes, elements' connectivity and boundary conditions


#FR = 20 # Frequency of the problem [Hz]
CW = 100 # Wave propagation speed [m/s]
#k = FR/CW # Wave number
fc = 0
finc = 0
FR1 = 22;
FR2 = 28;
nfreq = 10;
phi = complex(zeros(12,nfreq));
G = 0
for i = 1:1:nfreq
FR=FR1 + i*(FR2-FR1)/nfreq
println("Frequency =$(FR) rad/s.")
k = FR/CW;
# Gaussian quadrature - generation of points and weights [-1,1]
npg=6; # Number of integration points
qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
b1 = 1:nnos # Array containing all the indexes for nodes and elements which will be used for integration
println("Building A and b matrices using the traditional colocation BEM for constant elements.")
G,H = cal_GeH(NOS,NOS_GEO,ELEM,k,fc,qsi,w);
@time A,b = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])  # Builds A and B matrices using the collocation technique and applying the boundary conditions
#println("Tamanho de b = $(size(b))")
x = A\b # Solves the linear system
phi[:,i],qphi = monta_phieq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
println("Evaluating at domain points.")
@time phi_pint = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi[:,i],qphi,fc,finc,qsi,w,k) # Evaluates at internal (or external) points
end

