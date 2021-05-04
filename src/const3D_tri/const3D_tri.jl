# Boundary element method implementation for the Helmholtz equation using constant tridimensional triangular elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Module for the constant three-dimensional triangular element
# Contains the dependencies for the triangular element integration. The main function is const3D_tri.solve() which builds the influence matrices, applies the boundary conditions, solves the linear system and returns the value of the velocity potential and its flux at boundary and domain points.

module const3D_tri
using DelimitedFiles,  LinearAlgebra, Statistics,KrylovMethods
include("dep.jl") # Includes the dependencies
function solve(info,PONTOS_int,BCFace,k,verbose=false)
    NOS_GEO,ELEM,elemint,CDC = info
    NOS = mostra_geoTRI(NOS_GEO,ELEM) #Generate the physical nodes for constant elements
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    # If there's a defined BCFace, the boundary condition matrix is created. If no BCFace is defined, then the mesh reader has already built the boudnary conditions matrix.
    if isempty(BCFace) != true
	CDC = gera_CDC(ELEM,BCFace); #Monta a matriz de condicoes de contorno
    end
    # Gaussian quadrature - generation of points and weights [-1,1]
    npg=6; # Number of integration points
    qsi,w = Gauss_Legendre(0,1,npg) # Generation of the points and weights
    if verbose == true
        println("Building G and H matrices...")
	@time G,H,phi_inc = cal_GeH(NOS,NOS_GEO,ELEM,k,qsi,w,0) #Compute the G and H matrices
        println("Applying boundary conditions to build A and b for the linear system...")
        @time A,b = aplica_cdc(G,H,CDC) #Applies the boundary conditions and returns matrix A and vector b for the linear system
        println("Solving the linear system...")
        @time x = A\b # Solves the linear system
        println("Separating acoustic pressure from flux...")
        @time phi,q = monta_Teq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
        println("Solving for domain points.")
        @time phi_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
        @time dphidx,dphidy,dphidz=calc_q_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
    else
        G,H,phi_inc = cal_GeH(NOS,NOS_GEO,ELEM,k,qsi,w,0) #Compute the G and H matrices
        A,b = aplica_cdc(G,H,CDC) #Applies the boundary conditions and returns matrix A and vector b for the linear system
        x = A\b # Solves the linear system
        phi,q = monta_Teq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
        phi_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
        dphidx,dphidy,dphidz=calc_q_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
    end
    return phi,q,phi_pint,dphidz
end

function solveH(info,PONTOS_int,BCFace,k,verbose=false)
    NOS_GEO,ELEM,elemint,CDC = info
    NOS = mostra_geoTRI(NOS_GEO,ELEM) #Generate the physical nodes for constant elements
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    # If there's a defined BCFace, the boundary condition matrix is created. If no BCFace is defined, then the mesh reader has already built the boudnary conditions matrix.
    if isempty(BCFace) != true
	CDC = gera_CDC(ELEM,BCFace); #Monta a matriz de condicoes de contorno
    end
    Tree,block = cluster(NOS[:,2:4],floor(sqrt(length(NOS))),2)
    # Gaussian quadrature - generation of points and weights [-1,1]
    npg=12; # Number of integration points
    qsi,w = Gauss_Legendre(0,1,npg) # Generation of the points and weights
    A,b = Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,[0],qsi,w,CDC,k])
    x = gmres(vet->matvec(A,vet,block,Tree),b,5,tol=1e-5,maxIter=1000,out=0)
    phi,q = monta_Teq(CDC,x[1]) # Applies the boundary conditions to return the velocity potential and flux
    phi_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
    dphidx,dphidy,dphidz=calc_q_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
return phi,q,phi_pint,dphidz
end

function gera_vars(ELEM,CCFace,NOS_GEO)

    # Descri��o: Gera a matriz CDC a partir das condi��es de contorno das
    #            faces.
    # Autor:     Gustavo Gontijo
    #

    # Gera a matriz de CDC
    # CDC = [n�mero do elemento, tipo da CDC, valor da CDC no n� 1,...
    #                               valor da CDC no n� 2, valor da CDC no n� 3]
    nelem = size(ELEM,1)
    CDC = zeros(nelem,3)
    NOS=zeros(nelem,4)
    for i=1:nelem
        CDC[i,1] = i
        CDC[i,2] = CCFace[ELEM[i,5],2];
        CDC[i,3] = CCFace[ELEM[i,5],3];
        noselem = ELEM[i,2:4]
        X1=NOS_GEO[noselem[1],2]
        Y1=NOS_GEO[noselem[1],3]
        Z1=NOS_GEO[noselem[1],4]
        X2=NOS_GEO[noselem[2],2]
        Y2=NOS_GEO[noselem[2],3]
        Z2=NOS_GEO[noselem[2],4]
        X3=NOS_GEO[noselem[3],2]
        Y3=NOS_GEO[noselem[3],3]
        Z3=NOS_GEO[noselem[3],4]
        XM=(X1+X2+X3)/3;
        YM=(Y1+Y2+Y3)/3;
        ZM=(Z1+Z2+Z3)/3;
        NOS[i,1]=i
        NOS[i,2]=XM
        NOS[i,3]=YM
        NOS[i,4]=ZM
    end
    return CDC,NOS
end

end
