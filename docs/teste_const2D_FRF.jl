# Boundary element method implementation for the Helmholtz equation using constant bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com

using SpecialFunctions
using KrylovMethods
include("../src/const2D/dep.jl") # Includes the dependencies
include("../src/const2D/dad_1.jl") # Includes the data file containing the geometry and physical boundary conditions of the problem
include("../src/const2D/beminterp.jl") # H-Matrices using Lagrange polynomial interpolation
include("../src/const2D/ACA.jl") # H-Matrices using ACA
using PyPlot

fr1=2*pi*20*10; fr2=2*pi*20*100; npassos=2; passo = (fr2-fr1)/npassos;
t = zeros(npassos+1,4)
iter=1
for FR = fr1:passo:fr2
println("Frequency = $(FR/(2*pi)) Hz")
i =  50# Number of elements for half circle
#FR = 20*10 # Frequency of the problem [rad/s]
CW = 343 # Wave propagation speed [m/s]
k = FR/CW # Wave number
n_pint = 100;
PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc,phi_analytical = dad_1(i,n_pint,FR) # Geometric and physical information of the problem
println("$(i*4) Elements, simulating vibrating cylinder.")

# Gaussian quadrature - generation of points and weights [-1,1]
npg=6; # Number of integration points
qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
NOS_GEO,NOS,ELEM,CDC = format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg) # Apply the discretization technique and builds the problems matrices for the geometrical points, physical nodes, elements' connectivity and boundary conditions
nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
b1 = 1:nnos # Array containing all the indexes for nodes and elements which will be used for integration
println("Building A and b matrices using the traditional colocation BEM for constant elements.")
tic()
A,b = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])  # Builds A and B matrices using the collocation technique and applying the boundary conditions
t[iter,1] = toq()
tic()
x = A\b # Solves the linear system
t[iter,2] = toq()
phi,qphi = monta_phieq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
println("Evaluating values at internal points.")
@time phi_pint = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,finc,qsi,w,k) # Evaluates the value at internal (or external) points
println("Calculating the error.")
@time erro = abs((sum((phi_pint - phi_analytical).^2))/sum(phi_analytical))   # Calcula a norma em comparação com a solução analítica.
println("error = $erro %")
plot(PONTOS_int[:,2],real.(phi_analytical),label="A - $(round(FR/(2*pi),2)) Hz",marker="x")
plot(PONTOS_int[:,2],real.(phi_pint),label="BEM - $(round(FR/(2*pi),2)) Hz",marker="+",markersize=12)
#figure()
#surf(PONTOS_int[:,2],PONTOS_int[:,3],real.(phi_analytical),label="Analytical FR = $(FR)")
#figure()
#surf(PONTOS_int[:,2],PONTOS_int[:,3],real.(phi_pint),label="BEM FR = $(FR)")

## H-Matrix - Interpolation using Lagrange polynomial
println("Building Tree and blocks using H-Matrices.")
tic()
Tree,block = cluster(NOS[:,2:3],floor(sqrt(length(NOS))),2)
println("Building A and b matrices using H-Matrix with interpolation.")
Ai,bi = Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
t[iter,3] = toq()
tic()
xi = gmres(vet->matvec(Ai,vet,block,Tree),bi,5,tol=1e-5,maxIter=1000,out=0) #GMRES nas matrizes do ACA
t[iter,4] = toq()
iter+=1
phii,qphii = monta_phieq(CDC,xi[1]) # Applies the boundary conditions to return the velocity potential and flux
println("Evaluating values at internal points.")
@time phi_pinti = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phii,qphii,fc,finc,qsi,w,k) # Evaluates the value at internal (or external) points
println("Calculating the error.")
@time erro = abs((sum((phi_pinti - phi_analytical).^2))/sum(phi_analytical))   # Calcula a norma em comparação com a solução analítica.
println("error = $erro %")
@time erroi = abs((sum((phi_pinti - phi_pint).^2))/sum(phi_pint))   # Calcula a norma em comparação com a solução analítica.
println("error against conventional bie = $erroi %")
figure(1)
plot(PONTOS_int[:,2],real.(phi_pinti),label="H-BEM - $(round(FR/(2*pi),2)) Hz",marker="*")
end

legend(loc=0)
grid(1)
xlabel("distance [m]")
ylabel("particle velocity [m/s]")

figure()
title("Linear system setup")
plot(1:size(t,1),t[:,1],label="BEM")
plot(1:size(t,1),t[:,3],label="H-BEM")
ylabel("seconds")
grid(1)
legend()

#figure()
#title("Linear system solving")
#plot(1:size(t,1),t[:,2],label="BEM")
#plot(1:size(t,1),t[:,4],label="H-BEM")
#ylabel("seconds")
#grid(1)
#legend()
#figure(3)
#surf(PONTOS_int[:,2],PONTOS_int[:,3],real.(phi_pinti),label="H-BEM FR = $(FR)")
