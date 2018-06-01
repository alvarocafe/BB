# Boundary element method implementation for the Helmholtz equation using constant bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com

using SpecialFunctions
using KrylovMethods
include("dep.jl") # Includes the dependencies
include("dad_1.jl") # Includes the data file containing the geometry and physical boundary conditions of the problem
include("H_mat.jl") # H-Matrices support for building the cluster tree and blocks
include("beminterp.jl") # H-Matrices using Lagrange polynomial interpolation
include("ACA.jl") # H-Matrices using ACA

i =  10# Number of elements for half circle
FR = 20 # Frequency of the problem [Hz]
CW = 343 # Wave propagation speed [m/s]
k = FR/CW # Wave number
PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc,phi_analytical = dad_1(i,FR) # Geometric and physical information of the problem
println("$(i*4) Elements, simulating vibrating cylinder.")
# Gaussian quadrature - generation of points and weights [-1,1]
npg=8; # Number of integration points
qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
NOS_GEO,NOS,ELEM,CDC = format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg) # Apply the discretization technique and builds the problems matrices for the geometrical points, physical nodes, elements' connectivity and boundary conditions
nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
b1 = 1:nnos # Array containing all the indexes for nodes and elements which will be used for integration
println("Building A and b matrices using the traditional colocation BEM for constant elements.")
G,H=cal_GeH(NOS,NOS_GEO,ELEM,k,fc,qsi,w);
@time A,b = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])  # Builds A and B matrices using the collocation technique and applying the boundary conditions
#println("Tamanho de b = $(size(b))")
x = A\b # Solves the linear system
phi,qphi = monta_phieq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
println("Evaluating at domain points.")
@time phi_pint = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,finc,qsi,w,k) # Evaluates the value at internal (or external) points
println("Calculating the error.")
@time erro = abs((sum((phi_pint - phi_analytical).^2))/sum(phi_analytical))   # Calcula a norma em comparação com a solução analítica.
println("error = $erro %")

## H-Matrix - Interpolation using Lagrange polynomial
println("Building Tree and blocks using H-Matrices.")
@time Tree,block = cluster(NOS[:,2:3],floor(sqrt(length(NOS))),2)
println("Building A and b matrices using H-Matrix with interpolation.")
@time Ai,bi = Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
xi = gmres(vet->matvec(Ai,vet,block,Tree),bi,5,tol=1e-5,maxIter=1000,out=0) #GMRES nas matrizes do ACA
phii,qphii = monta_phieq(CDC,xi[1]) # Applies the boundary conditions to return the velocity potential and flux
println("Evaluating at domain points.")
@time phi_pinti = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phii,qphii,fc,finc,qsi,w,k) # Evaluates the value at internal (or external) points
println("Calculating the error.")
@time erro = abs((sum((phi_pinti - phi_analytical).^2))/sum(phi_analytical))   # Calcula a norma em comparação com a solução analítica.
println("error = $erro %")
@time erroi = abs((sum((phi_pinti - phi_pint).^2))/sum(phi_pint))   # Calcula a norma em comparação com a solução analítica.
println("error against conventional bie = $erroi %")
