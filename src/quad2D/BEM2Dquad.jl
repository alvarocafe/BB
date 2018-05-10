# Boundary element method implementation for the Helmholtz equation using quadratic bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Modules necessary: SpecialFunctions.jl, PyPlot.jl
workspace()
using SpecialFunctions
using PyPlot
include("dep.jl")	# Include the dependencies
#Input dos dados
include("dad.jl")	# Data entry functions
NOS, NOS, ELEM, CDC, PONTOS_int, phi_analytical, CW, FR, fc, finc = dad_quad_cyl(10)
#Geração dos pontos e pesos de Gauss
npg=6*6;
qsi,w = Gauss_Legendre(-1,1.00001,npg)
# Building matrices H and G
println("Building G and H matrices...")
G,H = cal_GeH(NOS,ELEM,CW,FR,fc,qsi,w) # Evaluates the influence matrices H and G
println("Soma da matriz G[1,:]: ",sum(G[1,:]))
println("Building A and b matrices for the linear system...")
A,b,T_PR = aplica_CDC(G,H,CDC,ELEM)
x = b\A

