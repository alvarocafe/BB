# Boundary element method implementation for the Helmholtz and Laplace
#equations using constant  bidimensional elements
# Authors: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com,
# Fernando Barreto Soares - fernando.bsoares3@gmail.com
# Contains the dependencies for the linear constant element integration.
#The main function is const2D.solve() which builds the influence matrices,
#applies the boundary conditions, solves the linear system and returns the
#value of the potential and its gradient at boundary and domain points.
using PyCall, Plots, JLD
include("../BEM_base.jl")

### Analytical solution
L = 10; # Square length
#n = 1; # First mode number (Note: modes i and j = 0, excitation only on the "z" axis)
phi_cup(k,x) = sin.(k.*x) #closed-open
phi_closed(k,x) = cos.(k.*x) #closed-closed BCs
q_closed(k,x) = k*sin.(k.*x)
k_closed(n,L) = (n/2)*2*pi/L; # Resonance wavenumber, k = ω/c, where c is the speed of propagation and ω is the angular frequency
k_res(n=1,L=1) = (2*n-1)*pi/(2*L) #k_cup

CW = 343; # Speed of sound in m/s
#k = k_res(n,L);
k1 = k_res(1,L); # Set the frequency in [rad/s]
k2 = k_res(2,L); # Set the frequency in [rad/s]
k3 = k_res(3,L); # Set the frequency in [rad/s]
#FR = k*CW;
inc = [0]; # There'll be no incident wave

# BEM modelling
# Gaussian quadrature - generation of points and weights [-1,1]
npg=6; # Number of integration points
qsi,w = const3D_tri.Gauss_Legendre(-1,1,npg) # Generation of the points and weights
# Python - mesh.io
meshio = pyimport("meshio")
# Build the domain points
L = 10; # Length of the cube
#n_pint = 40; # Number of domain points
n_pint = 5; # Number of domain points
PONTOS_int = zeros(n_pint,4);
delta = 0.1; # distance from both ends 
passo = (L-2*delta)/(n_pint-1);
for i = 1:n_pint
    PONTOS_int[i,:] = [i 0 0 delta+(i-1)*passo];
end
# Set the boundary conditions for each face. Faces 1 and 3 are perpendicular do "z" axis, 4 and 2 are perpendicular to "x" and 6 and 5 are perpendicular to "y".
#BCFace = [1. 1. 0.
#          2. 1. 0.
#          3. 0. 1.
#          4. 1. 0.
#          5. 1. 0.
#          6. 1. 0.];
#BCFace = [1. 0. -1.
#          2. 1. 0.
#          3. 0. 1.
#          4. 1. 0.
#          5. 1. 0.
#          6. 1. 0.];
#BCFace = [1. 0. 0.
#          2. 1. 0.
#          3. 0. 1.
#          4. 1. 0.
#          5. 1. 0.
#          6. 1. 0.];
#BCFace = [1. 1. -1.
#          2. 1. 0.
#          3. 1. 1.
#          4. 1. 0.
#          5. 1. 0.
#          6. 1. 0.];          
#BCFace = [1. 1. 0.
#          2. 1. 0.
#          3. 1. 1.
#          4. 1. 0.
#          5. 1. 0.
#          6. 1. 0.];
BCFace = [1. 0. 0.
          2. 1. 0.
          3. 1. 1.
          4. 1. 0.
          5. 1. 0.
          6. 1. 0.];
mshd = "./data/"
file = "cube_coarse.msh"
#t = [];
#Gnelem = [];
mesh = meshio.read(string(mshd,file))
NOS_GEO = [1:size(mesh.points,1) mesh.points]
nelem = size(mesh.cells_dict["triangle"],1)
ELEM = [1:nelem mesh.cells_dict["triangle"].+1 mesh.cell_data["gmsh:geometrical"][3]]
CDC,NOS = const3D_tri.gera_vars(ELEM,BCFace,NOS_GEO)
println("Malha ",file,", nelem = ",nelem)
elemint = []
info = [NOS_GEO,ELEM,elemint,CDC]

pontos_freq = 400; # Number os frequency points
#pontos_freq = 20; # Number os frequency points
kFRFF = zeros(pontos_freq);
#intervalo_freq = k; # Set the frequency range
intervalo_freq = (k3-k1)+1.3*k1
passo_freq = intervalo_freq/pontos_freq;
kFRFF[1]=(k2+k_closed(2,L))/2-intervalo_freq/2; # First frequency value.
#kFRFF[1]=k-intervalo_freq/2;
for i = 2:pontos_freq
  kFRFF[i] = kFRFF[i-1]+passo_freq;
end

TcFRF =[]
T_pintcFRF = zeros(ComplexF64,n_pint,length(kFRFF));
for i = 1:length(kFRFF)
  kFRF = kFRFF[i]
  Tc,qc,T_pintc,qzc = const3D_tri.solve(info,PONTOS_int,BCFace,kFRF,true)
  append!(TcFRF,Tc)
  #append!(T_pintcFRF,T_pintc)
  T_pintcFRF[:,i] = T_pintc;
end
#TccH,qccH,T_pintccH,qzccH = const3D_tri.solveH(info,PONTOS_int,BCFace,k,true)

save(string("./data/FRF00.11.jld"),"TcFRF",TcFRF,"T_pintcFRF",T_pintcFRF)

# Results plotting

omega=kFRFF*CW;

pyplot()

#use the scatter plot to obtain the FRF
#scatter(omega,real(T_pintcFRF[1,:]),xlabel=("Frequência [rad/s]"), ylabel=("u"), marker=("o"), fillstyle=("none"),label=("z = 0.10"))
#scatter(omega,real(T_pintcFRF[2,:]),xlabel=("Frequência [rad/s]"), ylabel=("u"), marker=("o"), fillstyle=("none"),label=("z = 2.55"))
#scatter(omega,real(T_pintcFRF[3,:]),xlabel=("Frequência [rad/s]"), ylabel=("u"), marker=("o"), fillstyle=("none"),label=("z = 5.00"))
#scatter(omega,real(T_pintcFRF[4,:]),xlabel=("Frequência [rad/s]"), ylabel=("u"), marker=("o"), fillstyle=("none"),label=("z = 7.45"))
#scatter(omega,real(T_pintcFRF[5,:]),xlabel=("Frequência [rad/s]"), ylabel=("u"), marker=("o"), fillstyle=("none"),label=("z = 9.90"))
#For absolute value
internal_point = 2; #select desired internal point position on the PONTOS_int array for absolute FRF analysis
absoluteT_pintcFRF=zeros(length(kFRFF));
for i = 1:length(kFRFF)
  absoluteT_pintcFRF[i] = abs(T_pintcFRF[internal_point,i]);
end
scatter(omega,absoluteT_pintcFRF,xlabel=("Frequência [rad/s]"), ylabel=("u"), marker=("o"), fillstyle=("none"),label=("z = 2.55"))
vline!([k_res(1,10)*CW, k_res(2,10)*CW, k_res(3,10)*CW],label=("k_res"))
#vline!([k_closed(1,10)*CW, k_closed(2,10)*CW, k_closed(3,10)*CW],label=("k_res"))

#use the regular plot to obtain the pressure distribution inside the cube at a given frequency
#plot(PONTOS_int[:,4],real(T_pintcFRF[:,390]/maximum(real(T_pintcFRF[:,390]))),xaxis=("z"),yaxis=("u"),title=("Solução numérica de u"),label=("MEC"))
#plot!(PONTOS_int[:,4],phi_closed(k_closed(3,10),PONTOS_int[:,4]),label=("Analítico"))
