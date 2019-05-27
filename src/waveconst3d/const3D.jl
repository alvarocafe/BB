using KrylovMethods, Statistics, LinearAlgebra, DelimitedFiles
include("dad.jl")
include("hmat.jl")
include("bem_functions.jl")
include("lermsh.jl")
#include("pre_proc.jl")
include("tree.jl")

filemsh = "VT_A.msh"
NOS_GEO,ELEM,elemint,CDC = lermsh(filemsh,3) #Read the mesh generated
# Build the domain points
L = 140; # Length of the vocal tract
n_pint = 40; # Number of domain points
#n_pint = 10; # Number of domain points
PONTOS_int = zeros(n_pint,4);
delta = 1; # distance from both ends 
passo = (L-2*delta)/(n_pint-1);
for i = 1:n_pint
    PONTOS_int[i,:] = [i 0 0 delta+(i-1)*passo];
end
# Set the boundary conditions for each face. Vowel /A/ model has 30 faces
BCFace = ones(30,3);
BCFace[:,3] .= 0;
BCFace[1,:] = [1 1 1]; # Neumann (flux = 1) to the Glotis
BCFace[30,:] = [30 0 0]; # Dirichlet (pressure = 0) to the mouth
# u,q,uint,qint = const3D_tri.solve(mshinfo,PONTOS_int,BCFace,k)

# H-BEM
# max_elem = Define máximo de nós em cada folha, tal que: max_elem/2 <= nós em cada folha < max_elem
# max_elem=floor(sqrt(2*length(NOS[:,1])))
max_elem=10
println("max_elem = $max_elem")
Tree,child,center_row,diam,inode,ileaf = cluster(NOS[:,2:3],max_elem)
ninterp=4 # Número de pontos de interpolação
#η =.4 # Coeficiente relacionado a admissibilidade dos blocos
η =.7 # Coeficiente relacionado a admissibilidade dos blocos
allow=checa_admiss(η,center_row,diam,inode,ileaf) # mostra quais blocos são admissíveis
block = blocks(Tree,child,allow) # Função que retorna os blocos admissiveis
hmati,bi =Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,qsi,w,CDC,FR,CW],ninterp)
xi = gmres(vet->matvec(hmati,vet,block,Tree),bi,5,tol=1e-8,maxIter=1000,out=0) 
T,q=monta_Teq(CDC,xi[1])
T_pint = calc_T_pint(PONTOS_int,NOS_GEO,ELEM,T,q,FR,CW,qsi,w,inc)
