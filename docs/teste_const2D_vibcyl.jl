# Boundary element method implementation for the Helmholtz equation using constant bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com

using SpecialFunctions
include("../src/const2D/dep.jl") # Includes the dependencies
include("../src/const2D/dad_1.jl") # Includes the data file containing the geometry and physical boundary conditions of the problem

ne = 20 # Number of elements for half cylinder
CW = 343 # Speed of sound
FR = 1000 # Frequency
k = FR/CW # Wavenumber
n_pint = 10 # 
raio = 0.1 # Radius of the vibrating cylinder
delta = 10 # maximum distance from cylinder
PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int,fc,finc,phi_analytical = dad_vibcyl(ne,k,n_pint,raio,delta) # Geometric and physical information of the problem

# Gaussian quadrature - generation of points and weights [-1,1]
npg=6; # Number of integration points
qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
NOS_GEO,NOS,ELEM,CDC = format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg) # Apply the discretization technique and builds the problems matrices for the geometrical points, physical nodes, elements' connectivity and boundary conditions
nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
b1 = 1:nnos # Array containing all the indexes for nodes and elements which will be used for integration
println("Building A and b matrices using the traditional colocation BEM for constant elements.")
A,b = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])  # Builds A and B matrices using the collocation technique and applying the boundary conditions
x = A\b # Solves the linear system
phi,qphi = monta_phieq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
println("Evaluating at domain points.")
@time phi_pint = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,finc,qsi,w,k) # Evaluates the value at internal (or external) points
println("Calculating the error.")
@time erro = abs((sum((phi_pint - phi_analytical).^2))/sum(phi_analytical))   # Calcula a norma em comparação com a solução analítica.
println("error = $erro %")

