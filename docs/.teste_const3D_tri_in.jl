using PyPlot
c = 600 # Speed of propagation in mm/s
FR = 1	# Frequency
k = FR/c;
#k=1
phi_closed(x,n,L,c) = cos.(n*pi*(x./L))
F_closed(n,L,c) = pi*n*c/L # Analytical resonance frequency in rad/s
include("../src/const3D_tri/dep.jl")
#NOS_GEO,ELEM,elemint,CDC = lermsh("dados/cylinder_in_tri_v.msh",3) #Read the mesh generated 
NOS_GEO = [1 0.0	 0.0	 0.0
 2 1.0	 0.0	 0.0
 3 0.0	 1.0	 0.0
 4 1.0	 1.0	 0.0
 5 0.0	 0.0	 1.0
 6 1.0	 0.0	 1.0
 7 0.0	 1.0	 1.0
 8 1.0	 1.0	 1.0]
ELEM = [1  1     4     2 1
2   1     3     4     1
3   1     6     5     2
4   1     2     6     2
5   2     8     6     3
6   2     4     8     3
7   3     8     4     4
8   3     7     8     4
9   1     7     3     5
10   1     5     7     5
11   5     8     7     6
12   5     6     8     6]
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
BCFace = [1 0 0
          2 1 0
          3 1 0
          4 1 0
	  5 1 0
    	  6 0 1] # This will act like a piston in the end of the cylinder so that the cavity is perturbed
CDC = gera_CDC(ELEM,BCFace)
# Now, the boundary conditions are applied to obtain the linear system
A,b = aplica_cdc(G,H,CDC) #Applies the boundary conditions and returns matrix A and vector b for the linear system
println("b = $(b)")
x = A\b # Solves the linear system
phi,q = monta_Teq(CDC,x); # Applies the boundary conditions to return the velocity potential and flux
# Now, the domain points will be created
n_pint = 100
L = 1;
PONTOS_int = zeros(n_pint,4)
for i = 1:n_pint
    PONTOS_int[i,:] = [i 0 0 (L/n_pint)*i]
end
# The velocity potential is obtained for domain points
phi_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
figure()
plot(PONTOS_int[:,4],real(phi_pint)./maximum(abs.(phi_pint[:,1])),label="BEM tri")
xlabel("distance [m]")
ylabel("normalized velocity potential")
grid(1)
plot(PONTOS_int[:,4],phi_closed(PONTOS_int[:,4],n,L,c),label="Mode n=$(n), F = $(F_closed(n,L,c)/(2*pi)) Hz")
legend()
