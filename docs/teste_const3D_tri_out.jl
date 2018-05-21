using PyPlot
close("all")
phi_closed(x,n,L,c) = cos.(n*pi*(x./L))
F_closed(n,L,c) = pi*n*c/L # Analytical resonance frequency in rad/s
c = 343*1000 # Speed of propagation in mm/s
n = 1 # Mode number
for n = 1:3
L = 100 # Length of the cylinder
FR = F_closed(n,L,c)	# Frequency
k = FR/c;
include("../src/const3D_tri/dep.jl")
NOS_GEO,ELEM,elemint,CDC = lermsh("../dados/cylinder_out_tri_v.msh",3) #Read the mesh generated 
NOS = mostra_geoTRI(NOS_GEO,ELEM);
#n = 1; L = 100; c = 343*1000;
#k = F_closed(n,L,c)/c
# Gaussian quadrature - generation of points and weights [-1,1]
npg=4; # Number of integration points
qsi,w = Gauss_Legendre(0,1,npg) # Generation of the points and weights
println("Building G and H matrices...")
@time G,H,phi_inc = cal_GeH(NOS,NOS_GEO,ELEM,k,qsi,w,0) #Compute the G and H matrices
# To apply the boundary conditions, when no BCFace is given to the lermsh command, all of the boundary conditions are assumed to be Neumann and null
# So, a new BCFace is built
BCFace = [1 1 0
          2 1 0
          3 1 0
          4 0 1] # This will act like a piston in the end of the cylinder so that the cavity is perturbed
CDC = gera_CDC(ELEM,BCFace)
# Now, the boundary conditions are applied to obtain the linear system
A,b = aplica_cdc(G,H,CDC) #Applies the boundary conditions and returns matrix A and vector b for the linear system
x = A\b # Solves the linear system
phi,q = monta_Teq(CDC,x); # Applies the boundary conditions to return the velocity potential and flux
# Now, the domain points will be created
using SpecialFunctions
raio = 5
n_pint = 100
PONTOS_int = zeros(n_pint,4)
dx = 10.1
dy = 10.1
dz = 50
passo = 1
phi_analytical = complex(zeros(size(PONTOS_int,1),1));
for i = 1:n_pint  # Para n_pint pontos internos
  PONTOS_int[i,:] = [i  dx+i*passo dy dz]
  bh1 = SpecialFunctions.besselh(1,2,(k)*raio);
  bh0 = SpecialFunctions.besselh(0,2,(k)*PONTOS_int[i,2]);
  phi_analytical[i,1] = (1/k).*(bh0./bh1);      #solucao anali•tica pela separacao de variaveis em coordenadas cilindricas da equacao de Helmholtz
end
# The velocity potential is obtained for domain points
phi_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
plot(PONTOS_int[:,2],real(phi_pint),label="BEM tri")
xlabel("distance [mm]")
ylabel("velocity potential")
grid(1)
plot(PONTOS_int[:,2],real(phi_analytical),label="Analytical solution FR = $(FR)")
end
legend()
