using PyPlot
close("all")
phi_closed(x,n,L,c) = cos.(n*pi*(x./L))
F_closed(n,L,c) = pi*n*c/L # Analytical resonance frequency in rad/s
c = 343*1000 # Speed of propagation in mm/s
inc = [0] # no incident waves
n = 1 # Mode number
#for n = 1:3 # Start of "n" for loop
L = 100 # Length of the cylinder
FR = F_closed(n,L,c)	# Frequency
k = FR/c;
include("../src/const3D_tri/dep.jl")
NOS_GEO,ELEM,elemint,CDC = lermsh("../dados/cylinder_in_tri_v.msh",3) #Read the mesh generated
NOS = mostra_geoTRI(NOS_GEO,ELEM);
#n = 1; L = 100; c = 343*1000;
#k = F_closed(n,L,c)/c
# Gaussian quadrature - generation of points and weights [-1,1]
npg=4; # Number of integration points
qsi,w = Gauss_Legendre(0,1,npg) # Generation of the points and weights
# To apply the boundary conditions, when no BCFace is given to the lermsh command, all of the boundary conditions are assumed to be Neumann and null
# So, a new BCFace is built
BCFace = [1 1 0
          2 1 0
          3 1 0
          4 0 1] # This will act like a piston in the end of the cylinder so that the cavity is perturbed
CDC = gera_CDC(ELEM,BCFace)
b1 = 1:size(NOS,1)
@time Ai,bi,Gi,Hi = cal_Aeb(b1,b1,[NOS,NOS_GEO,ELEM,k,qsi,w,inc,CDC])

@time G,H,phi_inc = cal_GeH(NOS,NOS_GEO,ELEM,k,qsi,w,inc) #Compute the G and H matrices
# Now, the boundary conditions are applied to obtain the linear system
A,b = aplica_cdc(G,H,CDC) #Applies the boundary conditions and returns matrix A and vector b for the linear system

println("Building G and H matrices...")
println("sum(sum(G-Gi)) = $(sum(sum(G-Gi)))")
println("sum(sum(H-Hi)) = $(sum(sum(H-Hi)))")
println("sum(sum(A-Ai)) = $(sum(sum(A-Ai)))")
println("sum(sum(b-bi)) = $(sum(sum(b-bi)))")
x = A\b # Solves the linear system
phi,q = monta_Teq(CDC,x); # Applies the boundary conditions to return the velocity potential and flux
# Now, the domain points will be created
n_pint = 100
L = 100;
PONTOS_int = zeros(n_pint,4)
for i = 1:n_pint
    PONTOS_int[i,:] = [i 0 0 (L/n_pint-0.01)*i]
end
# The velocity potential is obtained for domain points
phi_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
q_pint=calc_q_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
plot(PONTOS_int[:,4],real(phi_pint)./maximum(abs.(phi_pint[:,1])),label=L"$\phi$ BEM")
plot(PONTOS_int[:,4],real(q_pint)./maximum(abs.(q_pint[:,1])),label=L"$q$ BEM")
xlabel("distance [mm]")
ylabel("normalized velocity potential")
grid(1)
plot(PONTOS_int[:,4],phi_closed(PONTOS_int[:,4],n,L,c),label="Mode n=$(n), F = $(F_closed(n,L,c)/(2*pi)) Hz")
#end # End of "n" for loop
legend()
