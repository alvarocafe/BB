# Boundary element method implementation for the Laplace equation using constant tridimensional triangular elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Module for the constant three-dimensional triangular element
# Contains the dependencies for the triangular element integration. The main function is potconst3D.solve() which builds the influence matrices, applies the boundary conditions, solves the linear system and returns the value of the potential and its flux at boundary and domain points.

module potconst3D
include("hmat.jl")
include("bem_functions.jl")
include("tree.jl")
using SpecialFunctions, KrylovMethods, Statistics, LinearAlgebra

function solve(info,PONTOS_int,CCFace,k)
    NOS_GEO,ELEM,elemint,CDC = info
    nelem=size(ELEM,1)
    println("Nũmero de nós: $nelem")
    CDC,NOS = gera_vars(ELEM,CCFace,NOS_GEO);   # Gera a matriz de condições de contorno
    # Mostra a geometria do problema
    npg=4;
    qsi,w = Gauss_Legendre(-1,1,npg) # Gera os pontos e pesos de Gauss
    A,b = cal_Aeb([1:nelem],[1:nelem],[NOS,NOS_GEO,ELEM,xi,w,CDC,k])
    x = b\A
    T,q=monta_Teq(CDC,x)
    return T,q
end


function solveH(info,PONTOS_int,CCFace,k)
    
    NOS_GEO,ELEM,elemint,CDC = info
    nelem=size(ELEM,1)
    println("Nũmero de nós: $nelem")
    CDC,NOS = gera_vars(ELEM,CCFace,NOS_GEO);   # Gera a matriz de condições de contorno
    # Mostra a geometria do problema
    npg=4;
    qsi,w = Gauss_Legendre(-1,1,npg) # Gera os pontos e pesos de Gauss

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
    hmati,bi =Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,qsi,w,CDC,k],ninterp)

    xi = gmres(vet->matvec(hmati,vet,block,Tree),bi,5,tol=1e-8,maxIter=1000,out=0) 

    T,q=monta_Teq(CDC,xi[1])
    A=montacheia(hmati,block,Tree,nelem)
    println("A matriz foi dividida em $(length(block[:,3])) blocos")
    println("Dentre estes blocos, $(sum(block[:,3])) são aproximados por matrizes de baixo rank")
    println("Rank máximo: $(ninterp*ninterp)")
    println("norma = $(norm(T-NOS[:,4]))")

    return T,q
end
end
