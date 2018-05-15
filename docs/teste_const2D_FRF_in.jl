# Boundary element method implementation for the Helmholtz equation using constant bidimensional elements
# Test script for the closed acoustic duct example
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com

using SpecialFunctions
using KrylovMethods
include("../src/const2D/dep.jl") # Includes the dependencies
include("../src/const2D/dad_1.jl") # Includes the data file containing the geometry and physical boundary conditions of the problem
include("../src/const2D/beminterp.jl") # H-Matrices using Lagrange polynomial interpolation
include("../src/const2D/ACA.jl") # H-Matrices using ACA
using PyPlot


close("all")
c = 343*1000 # Speed of propagation in mm/s
F_closed(n,L,c) = pi*n*c/L # Analytical resonance frequency in rad/s
phi_closed(x,n,L,c) = cos.(n*pi*(x./L))	# Analytical mode shape
## Now, to define a new geometry, first one must declare the points, segments, boundary conditions, etc...
n = 1; # Mode number
L = 100; # Length of the duct in mm
d = 10; # Diameter of the duct in mm
c = 343*1000; # Speed of wave propagation in mm/s
FR = F_closed(n,L,c); # Resonance frequency
k = FR/c	# Wave number
points = [1 0 0; 2 L 0; 3 L d; 4 0 d];
segments = [1 1 2 0; 2 2 3 0; 3 3 4 0; 4 4 1 0];
ne = 10; # Number of element per segment
BCFace = [1 1 0; 2 1 0; 3 1 0; 4 0 1]; # Face 4 will act like a piston
fc = [0];
finc=[0];
# Now, the domain points will be created
n_pint = 100
PONTOS_int = zeros(n_pint,3)
for i = 1:n_pint
    PONTOS_int[i,:] = [i ((L-0.1)/n_pint)*i d/2]
end


#i1=2; i2=100; npassos=2; passo = floor((i2-i1)/npassos);
t = zeros(1,5)
iter=1
#for FR = fr1:passo:fr2
#for ne = i1:passo:i1+npassos*passo
ne = 50
MESH = [1 ne; 2 ne; 3 ne; 4 ne];
# Gaussian quadrature - generation of points and weights [-1,1]
npg=6; # Number of integration points
qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
NOS_GEO,NOS,ELEM,CDC = format_dad(points,segments,MESH,BCFace) # Apply the discretization technique and builds the problems matrices for the geometrical points, physical nodes, elements' connectivity and boundary conditions
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
@time erro = abs((sum((phi_pint - phi_closed(PONTOS_int[:,2],n,L,c)).^2))/sum(phi_closed(PONTOS_int[:,2],n,L,c)))   # Calcula a norma em comparação com a solução analítica.
println("error = $erro %")
plot(PONTOS_int[:,2],real.(phi_closed(PONTOS_int[:,2],n,L,c)),label="A - $(round(FR/(2*pi),2)) Hz",marker="x")
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
@time erro = abs((sum((phi_pinti - phi_closed(PONTOS_int[:,2],n,L,c)).^2))/sum(phi_closed(PONTOS_int[:,2],n,L,c)))   # Calcula a norma em comparação com a solução analítica.
println("error = $erro %")
@time erroi = abs((sum((phi_pinti - phi_pint).^2))/sum(phi_pint))   # Calcula a norma em comparação com a solução analítica.
println("error against conventional bie = $erroi %")
figure(1)
plot(PONTOS_int[:,2],real.(phi_pinti),label="H-BEM - $(round(FR/(2*pi),2)) Hz",marker="*")
#end

legend(loc=0)
grid(1)
xlabel("distance [mm]")
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
