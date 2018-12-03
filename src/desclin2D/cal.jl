function monta_GeH(ELEM,NOS,CDC,k,fc,qsi,w)
    # Evaluates matrices H, G and array q, which corresponds to the influences of the velocity potentials, flux and external sources. The discretized boundary integral equation is, in its matricial form, H phi = G qphi + q. 

    n_el = length(ELEM[:,1]);	# Number of elements
    n_nos = length(NOS[:,1]);	# Number of physical nodes

    # Allocating matrices H, G and array q
    H = complex(zeros(n_nos,n_nos));
    G = complex(zeros(n_nos,n_nos));
    q=complex(zeros(n_nos,1));
    #H = zeros(n_nos,n_nos)
    #G = zeros(n_nos,n_nos)
    #q = zeros(n_nos)
    # A = complex(zeros(n_nos,n_nos));
    # B = complex(zeros(n_nos,n_nos));
    for i = 1 : n_nos	# Loop over the source points (physical nodes)
        # Coordinates of the source point
        xd = NOS[i,2];	# x coordinate of the source point
        yd = NOS[i,3];	# y coordinate of the source point
        for j = 1 : n_el	# Loop over the elements
            # Numbering of the nodes of element j
            no1 = ELEM[j,2];	
            no2 = ELEM[j,3];	 
            # Coordinate of the nodes [x1,y1,x2,y2]
            x1 = NOS[no1,2];	y1 = NOS[no1,3];
            x2 = NOS[no2,2];	y2 = NOS[no2,3];

            if ((i == no1) || (i == no2))  # Evaluates the singular integration term when node i belongs to element j
                if(i==no1)	# In this case, the node is at the point described by qsi = -1/2, at a quarter of the way into the line
                    eta,Jt = telles(qsi,-1/2);
                    g,h =calc_gh_nsing(x1,y1,x2,y2,xd,yd,k,eta,w.*Jt); # Evaluates the singular terms using the kernel for non singular integration
                    #h =h;
                else	# In this case, the node is at the point described by qsi = 1/2, at a quarter of the way to the end of the line
                    eta,Jt = telles(qsi,1/2);
                    g,h =calc_gh_nsing(x1,y1,x2,y2,xd,yd,k,eta,w.*Jt); # Evaluates the singular terms using the kernel for non singular integration
                    #h =h;
                end
	    else
                # Non singular integration
                g,h = calc_gh_nsing(x1,y1,x2,y2,xd,yd,k,qsi,w); # Evaluates the non singular term
            end
	    # Allocates the terms of matrices H and G
            H[i,(2*j-1):2*j] = h;
            G[i,(2*j-1):2*j] = g;
            # Applies the boundary conditions
            #        if CDC[j,2]==0
            #            A[i,2*j-1:2*j] = -h;
            #            B[i,2*j-1:2*j] = -g;
            #        else
            #            A[i,2*j-1:2*j] = g;
            #            B[i,2*j-1:2*j] = h;
            #        end
        end
        if fc[1,1]!=0
            q[i]=calc_q(xd,yd,fc,k);
        else
            q[i]=0;
        end
    end

    # Evaluation of the singular terms of matrix H using the rigid body consideration.
    for m = 1 : n_nos # Loop over the nodes
        H[m,m] = 0; # Let the diagonal term be zero
        for n = 1 : n_nos	# Loop over the elements
            if n != m	# If the term is not on the diagonal
                H[m,m] = H[m,m] - H[m,n];	# Subtract the term from the diagonal
            end;
        end;
    end;

    return G,H,q
end
function Monta_Teq(x,CDC)

    # Function to reorder de vectors in accordance with boundary conditions
    # Separa fluxo e temperatura

    # ncdc = n�mero de linhas da matriz CDC
    # T = vetor que cont�m as temperaturas nos n�s
    # q = vetor que cont�m o fluxo nos n�s

    n_el=size(CDC,1);
    todos_valores=zeros(1,2*n_el);
    todos_tipos = zeros(1,2*n_el);
    for i=1:n_el
        todos_valores[2*i-1] = CDC[i,3];
        todos_valores[2*i] = CDC[i,5];
        todos_tipos[2*i-1] = CDC[i,2];
        todos_tipos[2*i] = CDC[i,4];
    end

    ncdc = 2*n_el;

    # T = complex(zeros(ncdc))
    # q = complex(zeros(ncdc))
    T = zeros(ncdc)
    q = zeros(ncdc)
    for i=1:ncdc # La�o sobre as condi��es de contorno
        tipoCDC=todos_tipos[i]; # Tipo da condi��o de contorno
        valorCDC=todos_valores[i]; # Valor da condi��o de contorno
        valorcalculado=x[i]; # Valor que antes era desconhecido
        if tipoCDC == 1 # Fluxo � conhecido
            T[i] = valorcalculado; # A temperatura � o valor calculado
            q[i] = valorCDC; # O fluxo � a condi�ao de contorno
        else # A temperatura � conhecida
            T[i] = valorCDC; # A temperatura � a condi�ao de contorno
            q[i] = valorcalculado; # O fluxo � o valor calculado
        end
    end


    return T,q
end
function calc_gh_nsing(x1,y1,x2,y2,xd,yd,k,qsi,w)
    # Evaluates the non singular terms of matrices G and H.
    npg=size(qsi,1); # Number of integration (Gaussian quadrature) points

    #h=[0. 0.]; # Initializes the matrix H terms
    #g=[0. 0.]; # Initializes the matrix G terms
    h = complex(zeros(2)) # Initializes the matrix H terms
    g = complex(zeros(2)) # Initializes the matrix G terms

    d=sqrt((x2-x1)^2+(y2-y1)^2); # Distance between the two points (half the length of the element)
    L = 2*d;	# Length of the element
    dgamadqsi=L/2; # Jacobian

    sx=(x2-x1)/d; # x component of the tangent vector
    sy=(y2-y1)/d; # y component of the tangent vector
    nx=sy; # x component x of the normal vector
    ny=-sx; # y component of the normal vector

    for kk=1:npg	# Loop over the integration points
        N1,N2=calc_fforma_d(qsi[kk]); # Evaluates the shape functions

        x=N1*x1+N2*x2; # Evaluates the x coordinate of the integration point
        y=N1*y1+N2*y2; # Evaluates the y coordinate of the integration point

        phiast,qast=calc_solfund(x,y,xd,yd,nx,ny,k); # Evaluate the fundamental solutions
        #phiast = 1;
        h=h+[N1; N2].*qast.*dgamadqsi.*w[kk]; # Integra��o da matriz h
        g=g+[N1; N2].*phiast.*dgamadqsi.*w[kk]; # Integra��o da matriz g
    end
    return g,h
end
function calc_g_sing(x1,y1,x2,y2,k)
    # Evaluates the singular term for matrix g (this subroutine is for the continuous linear element)
    L = sqrt((x1-x2)^2+(y1-y2)^2);
    g = L/(4*pi*k)*(3/2-log(L));
    return g
end

function aplica_CDC(G,H,CDC)
                                        # Applies the boundary conditions by exchanging the columns of matrices G and H

                                        n_el=size(CDC,1); # Number of boundary conditions
                                        # Now two new matrices will be built to apply the boundary conditions. The reason for this style of boundary condition collocation is that the boundary condition matrix stored is for the continuous linear element.
                                        todos_valores=zeros(1,2*n_el);	# Allocate matrix for all the values of the boundary conditions
                                        todos_tipos = zeros(1,2*n_el);	# Allocate matrix for all the types of boundary conditions
                                        for i=1:n_el	# Loop over the elements
                                            todos_valores[2*i-1] = CDC[i,3];
                                            todos_valores[2*i] = CDC[i,5];
                                            todos_tipos[2*i-1] = CDC[i,2];
                                            todos_tipos[2*i] = CDC[i,4];
                                        end

                                        ncdc = 2*n_el; # The number of boundary conditions is twice the number of elements, as there are two nodes for each element
                                        A=copy(H);
                                        B=copy(G);
                                        for i=1:ncdc # Loop over the boundary conditions
                                            tipoCDC = todos_tipos[i]; # Type of boundary condition
                                            if tipoCDC == 0 # The velocity potential is known
                                                colunaA=-A[:,i];
                                                A[:,i]=-B[:,i]; # Matrix A receives the column from matrix G
                                                B[:,i]=colunaA; # Matrix B receives the column from matrix H
                                            end
                                        end;

                                        b=B*todos_valores'; # vector b
                                        return A,b
                                    end
function  calc_solfund(x,y,xd,yd,nx,ny,k)
    # Evaluates the fundamental solution for the Laplace equation
    R=sqrt((x-xd)^2+(y-yd)^2); # Distance between the source and field points
    Rx=(x-xd); # x component of the distance
    Ry=(y-yd); # y component of the distance
    Tast=-1/(2*pi*k)*log(R); # Fundamental solution for the temperature
    qast=1/(2*pi)*(Rx*nx+Ry*ny)/R^2; # Fundamental solution for the flux
    return Tast,qast
end
# function  calc_solfund(x,y,xd,yd,nx,ny,k)
# Evaluates the fundamental solutions for the Helmholtz equation
#  r=sqrt((x-xd)^2+(y-yd)^2); # Raio (distância entre ponto fonte e ponto campo)
#  rx=(x-xd)/r; # Componente x do raio
#  ry=(y-yd)/r; # Componente y do raio
#  drdn=rx*nx+ry*ny;   #Componente do raio na direção normal
#  ZR=real(k*r);
#  Z=complex(0.,ZR);
#  F0C=SpecialFunctions.besselk(0,Z);
#  F1C=SpecialFunctions.besselk(1,Z);
#
#  qast=-(Z/r*drdn*F1C)/(2*pi); #Solução Fundamental da pressão acústica
#  Tast=F0C/(2*pi);    #Solução Fundamental do fluxo de pressão acústica
#  return Tast,qast
# end
