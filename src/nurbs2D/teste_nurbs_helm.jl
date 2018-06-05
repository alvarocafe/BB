# Boundary element method implementation for the Laplace equation using NURBS bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Necessary Modules: SpecialFunctions.jl
using SpecialFunctions
using PyPlot
include("dep.jl") # Includes the dependencies
include("dad_1.jl") # Includes the data file containing the geometry and physical boundary conditions of the problem
# Characteristics of the problem: Square domain with imposed temperature in two opposite faces and imposed null temperature flux at the other two faces. 
k=0	# wavenumber
collocCoord,nnos,crv,dcrv,CDC,E = dad_helm()# Geometric and physical information of the problem

#Building the problems matrices
H, G = calcula_iso(collocCoord,nnos,crv,dcrv,E,k) # Influence matrices
A,b= aplica_CDCiso(G,H,CDC,E);	# Applying boundary conditions
x=A\b; # Evaluating unknown values
Tc,qc=monta_Teqiso(CDC,x); # Separating temperature from flux
# Applying NURBS basis functions to the values of temperature and flux
T=E*Tc;
q=E*qc;

# Domain points
PONTOS_int = [1 3 3]
fc = 0; finc = 0;
Hp,Gp,phi_pint = calc_phi_pint_nurbs(PONTOS_int,collocCoord,nnos,crv,dcrv,k,Tc,qc);

# Graphics
close("all")
mostra_geo(crv)
plot(collocCoord[:,1],collocCoord[:,2],marker="s",markersize=10,linestyle="none",color="blue",label = "Physical points (Nodes)")
axis("equal")
grid(1)
PyPlot.xlabel("x",fontsize="12.0")
PyPlot.ylabel("y",fontsize="12.0")
title("NURBS model",fontsize="16.0")
legend(fontsize="14.0",loc="best")
# Plot domain point solution and analytical solution
#figure()
#plot(PONTOS_int[:,2],real(phi_analytical),label="Analytical")
#plot(PONTOS_int[:,2],real(phi_pint),label="IGA BEM")
#legend()
