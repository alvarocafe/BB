# Boundary element method implementation for the Helmholtz equation using NURBS bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Necessary Modules: SpecialFunctions.jl
using SpecialFunctions
using KrylovMethods
using PyPlot
include("dep.jl") # Includes the dependencies
include("dad_1.jl") # Includes the data file containing the geometry and physical boundary conditions of the problem
# Characteristics of the problem: vibrating cylinder in an infinite acoustic medium.
# The geometry of the problem is a circle of diameter 1 vibrating in an infinite acoustic medium and the pressure field is defined by the Helmholtz equation. A concentrated acoustic source and incident plane waves can be added.

i =  1# Number of elements for half circle
FR = 20 # Frequency of the problem [Hz]
CW = 343 # Wave propagation speed [m/s]
k = FR/CW # Wave number
collocCoord,nnos,crv,dcrv,CDC,E = dad_iso(10)# Geometric and physical information of the problem
#println("$(i*4) Elements, simulating vibrating cylinder.")
# Gaussian quadrature - generation of points and weights [-1,1]
npg=16; # Number of integration points
qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
H, G = calcula_iso(collocCoord,nnos,crv,dcrv,1)
H=H+E/2;
println("Para o calculo do perimetro, pi = $(sum(G[1,:])*2)")
#println(sum(H[9,:])) # Mostra a soma da linha para o elemento de comparação
A,b= aplica_CDCiso(G,H,CDC,E);
x=A\b; # Calcula o vetor x
Tc,qc=monta_Teqiso(CDC,x); # Separa temperatura e fluxo
T=E*Tc;
q=E*qc;
close("all")
mostra_geo(crv)
plot(collocCoord[:,1],collocCoord[:,2],marker="s",markersize=10,linestyle="none",color="blue",label = "Physical points (Nodes)")
axis("equal")
grid(1)
PyPlot.xlabel("x",fontsize="12.0")
PyPlot.ylabel("y",fontsize="12.0")
title("NURBS element",fontsize="16.0")
legend(fontsize="14.0",loc="best")
PyPlot.xlim(-0.1, 1.1)
PyPlot.ylim(-0.6, 0.1)
#savefig("docs/discretização/figuras/nurbs1.pdf",transparent = true)
