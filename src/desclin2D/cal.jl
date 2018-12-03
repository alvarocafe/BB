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

    # ncdc = nï¿½mero de linhas da matriz CDC
    # T = vetor que contï¿½m as temperaturas nos nï¿½s
    # q = vetor que contï¿½m o fluxo nos nï¿½s

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
    for i=1:ncdc # Laï¿½o sobre as condiï¿½ï¿½es de contorno
        tipoCDC=todos_tipos[i]; # Tipo da condiï¿½ï¿½o de contorno
        valorCDC=todos_valores[i]; # Valor da condiï¿½ï¿½o de contorno
        valorcalculado=x[i]; # Valor que antes era desconhecido
        if tipoCDC == 1 # Fluxo ï¿½ conhecido
            T[i] = valorcalculado; # A temperatura ï¿½ o valor calculado
            q[i] = valorCDC; # O fluxo ï¿½ a condiï¿½ao de contorno
        else # A temperatura ï¿½ conhecida
            T[i] = valorCDC; # A temperatura ï¿½ a condiï¿½ao de contorno
            q[i] = valorcalculado; # O fluxo ï¿½ o valor calculado
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
        h=h+[N1; N2].*qast.*dgamadqsi.*w[kk]; # Integraï¿½ï¿½o da matriz h
        g=g+[N1; N2].*phiast.*dgamadqsi.*w[kk]; # Integraï¿½ï¿½o da matriz g
    end
    return g,h
end
function calc_g_sing(x1,y1,x2,y2,k)
    # Evaluates the singular term for matrix g (this subroutine is for the continuous linear element)
    L = sqrt((x1-x2)^2+(y1-y2)^2);
    g = L/(4*pi*k)*(3/2-log(L));
    return g
    function T_pint=calc_T_pint(ELEM,NOS,PONTOS_INT,T,q,k,fc)
        % Calcula a temperatura nos pontos internos
        n_pint=length(PONTOS_INT(:,1)); % Numero de pontos internos
        n_elem=length(T); % Numero de elementos
        H_int=zeros(n_pint,n_elem);
        G_int=zeros(n_pint,2*n_elem);
        q_fc=zeros(n_pint,1);

        npg=16;
        [qsi,w]=Gauss_Legendre(-1,1,npg);


        for i=1:n_pint % Laço sobre os pontos internos
            x_fonte=PONTOS_INT(i,2); % Coordenada x do ponto fonte
                y_fonte=PONTOS_INT(i,3); % Coordenada y do ponto fonte
                    for j=1:n_elem
                        x1=NOS(ELEM(j,2),2); % Coordenada x do inicio do elemento
                            y1=NOS(ELEM(j,2),3); % Coordenada y do inicio do elemento
                                x2=NOS(ELEM(j,3),2); % Coordenada x do final do elemento
                                    y2=NOS(ELEM(j,3),3); % Coordenada y do final do elemento
                                        [g,h]=calc_gh_nsing(x1,y1,x2,y2,x_fonte,y_fonte,k,qsi,w); % Chama a
                                        % functio para cálculo de H e G quando o
                                        % ponto fonte nao pertence ao elemento
                                        for nolocal = 1 : 2
                                            noglobal = ELEM(j,nolocal+1); %Índice da matriz global H
                                            H_int(i,noglobal) = H_int(i,noglobal) + h(nolocal);
                                        end;
                                        G_int(i,2*j-1:2*j) = g;
                                    end
                                        q_fc(i)=calc_q(x_fonte,y_fonte,fc,k);
                                    end

                                    T_pint=H_int*T-G_int*q'+q_fc; % Vetor que contem a temperatura nos pontos internosfunction calc_q(xd,yd,fc,k)
                                    # Evaluates the influence from concentrated sources
                                    valor_fonte=fc[1];
                                    x_fonte=fc[2];
                                    y_fonte=fc[3];
                                    R=sqrt((x_fonte-xd)^2+(y_fonte-yd)^2); # Distance between the source and field poitns
                                    Tast=-1/(2*pi*k)*log(R); # Fundamental solution for the temperature
                                    q=valor_fonte*Tast;
                                    return q
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
                                    function [Dx,Dy,Sx,Sy]=calc_DeS(x1,y1,x2,y2,xd,yd,k,qsi,w)

                                        npg=length(qsi); % Número de pontos de Gauss

                                        Dx=[0 0]; % Inicializa a matriz h do elemento
                                            Dy=[0 0]; % Inicializa a matriz g do elemento
                                                Sx=[0 0]; % Inicializa a matriz h do elemento
                                                    Sy=[0 0]; % Inicializa a matriz g do elemento

                                                        L=sqrt((x2-x1)^2+(y2-y1)^2); % Comprimento do elemento
                                                            dgamadqsi=L/2; % Jacobiano

                                                            sx=(x2-x1)/L; % Componente x do vetor tangente
                                                                sy=(y2-y1)/L; % Componente y do vetor tangente
                                                                    nx=sy; % Componente x do vetor normal
                                                                        ny=-sx; % Componente y do vetor normal

                                                                            for kk=1:npg
                                                                                [N1,N2]=calc_fforma(qsi(kk)); % Calcula as funções de forma

                                                                                x=N1*x1+N2*x2; % Calcula a coordenada x do ponto de integração
                                                                                    y=N1*y1+N2*y2; % Calcula a coordenada y do ponto de integração
                                                                                        
                                                                                        [dTdx,dTdy,dqdx,dqdy]=calc_dsolfund(x,y,xd,yd,nx,ny,k); % Calcula 
                                                                                        % as derivadas das soluções fundamentais
                                                                                        Dx=Dx+dTdx*[N1,N2]*dgamadqsi*w(kk); % Integral da matriz H
                                                                                        Dy=Dy+dTdy*[N1,N2]*dgamadqsi*w(kk); % Integral da matriz H
                                                                                        Sx=Sx+dqdx*[N1,N2]*dgamadqsi*w(kk); % Integral da matriz G
                                                                                        Sy=Sy+dqdy*[N1,N2]*dgamadqsi*w(kk); % Integral da matriz G
                                                                                        
                                                                                        endfunction [dqdx,dqdy]=calc_dq(xd,yd,fc,k)

                                                                                        valor_fonte=fc(1);
                                                                                        x_fonte=fc(2);
                                                                                        y_fonte=fc(3);

                                                                                        %Calcula as soluções fundamentais


                                                                                        rx=(x_fonte-xd); % Componente x do raio
                                                                                            ry=(y_fonte-yd); % Componente y do raio
                                                                                                r=sqrt(rx^2+ry^2);
                                                                                                dTdx=1/(2*pi*k*r^2)*rx; % Solução fundamental da temperatura
                                                                                                dTdy=1/(2*pi*k*r^2)*ry; % Solução fundamental da temperatura

                                                                                                dqdx=valor_fonte*dTdx;
                                                                                                dqdy=valor_fonte*dTdy;function [dTdx,dTdy]=calc_dT_pint(ELEM,NOS,PONTOS_INT,T,q,k,fc)
                                                                                                    % Calcula a temperatura nos pontos internos
                                                                                                    n_pint=length(PONTOS_INT(:,1)); % Numero de pontos internos
                                                                                                    n_elem=length(T); % Numero de elementos
                                                                                                    Sx_int=zeros(n_pint,n_elem);
                                                                                                    Sy_int=zeros(n_pint,n_elem);
                                                                                                    Dx_int=zeros(n_pint,2*n_elem);
                                                                                                    Dy_int=zeros(n_pint,2*n_elem);
                                                                                                    dqdx_fc=zeros(n_pint,1);
                                                                                                    dqdy_fc=zeros(n_pint,1);

                                                                                                    npg=16;
                                                                                                    [qsi,w]=Gauss_Legendre(-1,1,npg);


                                                                                                    for i=1:n_pint % Laço sobre os pontos internos
                                                                                                        x_fonte=PONTOS_INT(i,2); % Coordenada x do ponto fonte
                                                                                                            y_fonte=PONTOS_INT(i,3); % Coordenada y do ponto fonte
                                                                                                                for j=1:n_elem
                                                                                                                    x1=NOS(ELEM(j,2),2); % Coordenada x do inicio do elemento
                                                                                                                        y1=NOS(ELEM(j,2),3); % Coordenada y do inicio do elemento
                                                                                                                            x2=NOS(ELEM(j,3),2); % Coordenada x do final do elemento
                                                                                                                                y2=NOS(ELEM(j,3),3); % Coordenada y do final do elemento
                                                                                                                                    
                                                                                                                                    [dx,dy,sx,sy]=calc_DeS(x1,y1,x2,y2,x_fonte,y_fonte,k,qsi,w); % Chama a
                                                                                                                                    % functio para cálculo de H e G quando o
                                                                                                                                    % ponto fonte nao pertence ao elemento
                                                                                                                                    for nolocal = 1 : 2
                                                                                                                                        noglobal = ELEM(j,nolocal+1); %Índice da matriz global H
                                                                                                                                        Sx_int(i,noglobal) = Sx_int(i,noglobal) + sx(nolocal);
                                                                                                                                        Sy_int(i,noglobal) = Sy_int(i,noglobal) + sy(nolocal);
                                                                                                                                    end;
                                                                                                                                    Dx_int(i,2*j-1:2*j) = dx;
                                                                                                                                    Dy_int(i,2*j-1:2*j) = dy;
                                                                                                                                end
                                                                                                                                    [dqdx_fc(i),dqdy_fc(i)]=calc_dq(x_fonte,y_fonte,fc,k);    
                                                                                                                                end

                                                                                                                                dTdx=Dx_int*q'-Sx_int*T-dqdx_fc; % Vetor que contem a temperatura nos pontos internos
                                                                                                                                dTdy=Dy_int*q'-Sy_int*T-dqdy_fc; % Vetor que contem a temperatura nos pontos internosfunction T_pint=calc_T_pint(ELEM,NOS,PONTOS_INT,T,q,k,fc)
                                                                                                                                % Calcula a temperatura nos pontos internos
                                                                                                                                n_pint=length(PONTOS_INT(:,1)); % Numero de pontos internos
                                                                                                                                n_elem=length(T); % Numero de elementos
                                                                                                                                H_int=zeros(n_pint,n_elem);
                                                                                                                                G_int=zeros(n_pint,2*n_elem);
                                                                                                                                q_fc=zeros(n_pint,1);

                                                                                                                                npg=16;
                                                                                                                                [qsi,w]=Gauss_Legendre(-1,1,npg);


                                                                                                                                for i=1:n_pint % Laço sobre os pontos internos
                                                                                                                                    x_fonte=PONTOS_INT(i,2); % Coordenada x do ponto fonte
                                                                                                                                        y_fonte=PONTOS_INT(i,3); % Coordenada y do ponto fonte
                                                                                                                                            for j=1:n_elem
                                                                                                                                                x1=NOS(ELEM(j,2),2); % Coordenada x do inicio do elemento
                                                                                                                                                    y1=NOS(ELEM(j,2),3); % Coordenada y do inicio do elemento
                                                                                                                                                        x2=NOS(ELEM(j,3),2); % Coordenada x do final do elemento
                                                                                                                                                            y2=NOS(ELEM(j,3),3); % Coordenada y do final do elemento
                                                                                                                                                                [g,h]=calc_gh_nsing(x1,y1,x2,y2,x_fonte,y_fonte,k,qsi,w); % Chama a
                                                                                                                                                                % functio para cálculo de H e G quando o
                                                                                                                                                                % ponto fonte nao pertence ao elemento
                                                                                                                                                                for nolocal = 1 : 2
                                                                                                                                                                    noglobal = ELEM(j,nolocal+1); %Índice da matriz global H
                                                                                                                                                                    H_int(i,noglobal) = H_int(i,noglobal) + h(nolocal);
                                                                                                                                                                end;
                                                                                                                                                                G_int(i,2*j-1:2*j) = g;
                                                                                                                                                            end
                                                                                                                                                                q_fc(i)=calc_q(x_fonte,y_fonte,fc,k);
                                                                                                                                                            end

                                                                                                                                                            T_pint=H_int*T-G_int*q'+q_fc; % Vetor que contem a temperatura nos pontos internos
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
#  r=sqrt((x-xd)^2+(y-yd)^2); # Raio (distÃ¢ncia entre ponto fonte e ponto campo)
#  rx=(x-xd)/r; # Componente x do raio
#  ry=(y-yd)/r; # Componente y do raio
#  drdn=rx*nx+ry*ny;   #Componente do raio na direÃ§Ã£o normal
#  ZR=real(k*r);
#  Z=complex(0.,ZR);
#  F0C=SpecialFunctions.besselk(0,Z);
#  F1C=SpecialFunctions.besselk(1,Z);
#
#  qast=-(Z/r*drdn*F1C)/(2*pi); #SoluÃ§Ã£o Fundamental da pressÃ£o acÃºstica
#  Tast=F0C/(2*pi);    #SoluÃ§Ã£o Fundamental do fluxo de pressÃ£o acÃºstica
#  return Tast,qast
# end
function  [dTdx,dTdy,dqdx,dqdy]=calc_dsolfund(x,y,xd,yd,nx,ny,k) 
    %Calcula as soluções fundamentais

    r=sqrt((x-xd)^2+(y-yd)^2); % Raio (distância entre ponto fonte e 
                                       % ponto campo)
    rx=(x-xd); % Componente x do raio
        ry=(y-yd); % Componente y do raio

            dTdx=1/(2*pi*k*r^2)*rx; % Solução fundamental da temperatura
            dTdy=1/(2*pi*k*r^2)*ry; % Solução fundamental da temperatura

            dqdx=(nx*(rx^2 - ry^2) + 2*ny*rx*ry)/(2*pi*r^4);
            dqdy=(ny*(-rx^2 +ry^2) + 2*nx*rx*ry)/(2*pi*r^4);
