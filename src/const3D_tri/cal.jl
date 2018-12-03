function cal_Aeb(b1,b2,arg)
    NOS, NOS_GEO, ELEM, k, CDC,qsi,w,qsi_tri,w_tri = arg;
    nnos = length(b1);
    nelem = length(b2)
    G=zeros(nnos,nelem); 
    H=zeros(nnos,nelem);     
    A = complex(zeros(nnos,nelem));
    b = complex(zeros(nnos, 1));
    ci=0
    for i in b1
        ci+=1
        xd=NOS[i,2]; 
        yd=NOS[i,3]; 
        zd=NOS[i,4];
        
        cj=0
        for j in b2
            cj+=1
	    tipoCDC = CDC[j,2]; 
	    
	    valorCDC = CDC[j,3];
	    nos = ELEM[j,2:4];		

            no1=ELEM[j,2]; 
            no2=ELEM[j,3]; 
            no3=ELEM[j,4]; 

            x1=NOS_GEO[no1,2]; 
            y1=NOS_GEO[no1,3]; 
            z1=NOS_GEO[no1,4]; 

            x2=NOS_GEO[no2,2]; 
            y2=NOS_GEO[no2,3]; 
            z2=NOS_GEO[no2,4]; 

            x3=NOS_GEO[no3,2]; 
            y3=NOS_GEO[no3,3]; 
            z3=NOS_GEO[no3,4]; 

            n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); 
            
            if i==j 
                g,h=calcula_HeGs_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsi,w,k); 
            else 
                g,h=calcula_HeGns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi_tri,w_tri,k); 
            end            
            if tipoCDC == 0 
        	A[ci,cj] = -g; 	
                if valorCDC == 0
                    b[ci,1] = b[ci,1] - 0; 
                else
                    b[ci,1] = b[ci,1] - h*valorCDC; 
                end                
            else
                A[ci,cj] = +h; 	
                if valorCDC == 0
                    b[ci,1] = b[ci,1] + 0; 
                else
                    b[ci,1] = b[ci,1] + g*valorCDC; 
                end
            end
        end
    end
    return A,b
end
function aplica_cdc(G,H,CDC)
# Aplica as condi��es de contorno trocando as colunas das matrizes H e G
ncdc = length(CDC[:,1]); # n�mero de linhas da matriz CDC
A=copy(H);
B=copy(G);
for i=1:ncdc # La�o sobre as condi��es de contorno
    tipoCDC = CDC[i,2]; # Tipo da condi��o de contorno
    if tipoCDC == 0 # A temperatura � conhecida
        colunaA=-A[:,i]; # Coluna da matriz H que ser� trocada
        A[:,i]=-B[:,i]; # A matriz H recebe a coluna da matriz G
        B[:,i]=colunaA; # A mstriz G recebe a coluna da matriz H
    end
end;

#valoresconhecidos=complex(CDC[:,3],CDC[:,4]); # Valores das condicoes de contorno
valoresconhecidos = CDC[:,3];
b=B*valoresconhecidos; # vetor b

return A,b
end

function telles(gamm,eet)

eest = eet^2 - 1;
term1 = eet*eest + abs(eest);
if term1 < 0
    term1 = (-term1)^(1/3);
    term1 = -term1;
else
    term1 = term1^(1/3);
end

term2 = eet*eest - abs(eest);
if term2 < 0
    term2 = (-term2)^(1/3);
    term2 = -term2;
else
    term2 = term2^(1/3);
end
GAMM = term1 + term2 + eet;


Q = 1 + 3*GAMM^2;
A = 1/Q;
B = -3*GAMM/Q;
C = 3*GAMM^2/Q;
D = -B;

eta = A*gamm.^3 + B*gamm.^2 + C*gamm + D;
Jt = 3*A*gamm.^2 + 2*B*gamm + C;
return  eta,Jt
end
function calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k)
#integra��o n�o singular
n_pint=length(qsi); # N�mero de pontos de integra��o.
g=complex(0,0); # Inicializa o somatorio de g
h=complex(0,0); # Inicializa o somatorio de h
N = zeros(3)
for l=1:n_pint # La�o sobre os pontos de integra��o
    for m=1:n_pint # La�o sobre os pontos de integra��o
	xi = (1 - qsi[l])*qsi[m]
        N =calc_fformatri(xi,qsi[l]); #  fun��es de forma
	#dN = calc_dfforma(xi,qsi[m])
        x=N[1]*x1+N[2]*x2+N[3]*x3; # coordenada x do ponto de integra��o
        y=N[1]*y1+N[2]*y2+N[3]*y3 # coordenada y do ponto de integra��o
        z=N[1]*z1+N[2]*z2+N[3]*z3; # coordenada z do ponto de integra��o
	J = calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,xi,qsi[m]);# jacobiano
        #Tast,qast=calc_solfund(x,y,z,xd,yd,zd,n,FR,CW); # Solu��es
        #  fundamentais
	Tast,qast=calc_solfund(x,y,z,xd,yd,zd,n,k); # Solu��es
        #  fundamentais
	#Tast = 1 #To test the numerical integration, the value of the fundamental solutions is set to 1. In this special case, the value of the sum of the rows of matrices H or G should equal the surface area of the model.
	#qast = 1        
	g=g+Tast*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz G
        h=h+qast*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz H
    end
end
return g,h
end

function calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsi,w,k)

# Integra��o singular das matrizes H e G. No caso da matriz G, o elemento
# triangular � decomposto em tr�s quadril�teros degenerados na forma
# de tri�ngulos. Os dois primeiros n�s destes quadril�teros s�o
# coincidentes (formam um s� v�rtice do tri�ngulo) e correspondem ao ponto
# onde existe a singularidade, ou seja, ao ponto fonte que se localiza no
# centr�ide do elemento. Isto faz com que haja uma concentra��o de pontos
# de integra��o junto � singularidade, al�m do jacobiano ser igual a zero
# na singularidade. No caso da matriz H, a integra��o � anal�tica e sempre
# ser� igual a 1/2.

npg=length(qsi); # N�mero de pontos de integra��o
g = 0.0; # inicializa��o da matriz G
Tast = complex(0,0);
qast = complex(0,0);
for kk=1:3
    x1t=xd; # coordenada x do primeiro n� do quadrilatero desgenerado
    y1t=yd; # coordenada y do primeiro n� do quadrilatero desgenerado
    z1t=zd; # coordenada z do primeiro n� do quadrilatero desgenerado
    x2t=xd; # coordenada x do segundo n� do quadrilatero desgenerado
    y2t=yd; # coordenada y do segundo n� do quadrilatero desgenerado
    z2t=zd; # coordenada z do segundo n� do quadrilatero desgenerado

    if(kk==1) # Terceiro e quarto n�s do primeiro quadril�tero desgenerado
        x3t=x1; # coordenada x do terceiro n� do quadrilatero desgenerado
        y3t=y1; # coordenada y do terceiro n� do quadrilatero desgenerado
        z3t=z1; # coordenada z do terceiro n� do quadrilatero desgenerado
        x4t=x2; # coordenada x do quarto n� do quadrilatero desgenerado
        y4t=y2; # coordenada y do quarto n� do quadrilatero desgenerado
        z4t=z2; # coordenada z do quarto n� do quadrilatero desgenerado

    elseif(kk==2) # Terceiro e quarto n�s do segundo quadril�tero
                                  # desgenerado
        x3t=x2; # coordenada x do terceiro n� do quadrilatero desgenerado
        y3t=y2; # coordenada y do terceiro n� do quadrilatero desgenerado
        z3t=z2; # coordenada z do terceiro n� do quadrilatero desgenerado
        x4t=x3; # coordenada x do quarto n� do quadrilatero desgenerado
        y4t=y3; # coordenada y do quarto n� do quadrilatero desgenerado
        z4t=z3; # coordenada z do quarto n� do quadrilatero desgenerado

    elseif(kk==3) # Terceiro e quarto n�s do terceiro quadril�tero
                                  # desgenerado
        x3t=x3; # coordenada x do terceiro n� do quadrilatero desgenerado
        y3t=y3; # coordenada y do terceiro n� do quadrilatero desgenerado
        z3t=z3; # coordenada z do terceiro n� do quadrilatero desgenerado
        x4t=x1; # coordenada x do quarto n� do quadrilatero desgenerado
        y4t=y1; # coordenada y do quarto n� do quadrilatero desgenerado
        z4t=z1; # coordenada z do quarto n� do quadrilatero desgenerado

    end

    for ii = 1: npg; # la�o sobre a primeira vari�vel de integra��o
        for jj = 1: npg # la�o sobre a segunda vari�vel de integra��o
            N1,N2,N3,N4 = calc_fforma(qsi[ii],qsi[jj]); # Fun��es de
                                                       #  forma
            xc = x1t*N1+x2t*N2+x3t*N3+x4t*N4; # coordenada x do
            # ponto campo
            yc = y1t*N1+y2t*N2+y3t*N3+y4t*N4; # coordenada y do
            # ponto campo
            zc = z1t*N1+z2t*N2+z3t*N3+z4t*N4; # coordenada z do
            # ponto campo
            J = calc_jacobiano_quad(x1t,y1t,z1t,x2t,y2t,z2t,x3t,y3t,z3t,x4t,y4t,z4t,qsi[ii],qsi[jj]); # jacobiano(varia ao longo
                         #  do elemento desgenerado)
            Tast,qast = calc_solfund(xd,yd,zd, xc, yc, zc, [0 0 0], k); # Sol.
            # fudamental de pressão acústica
            g = g + w[jj]*w[ii]*J*Tast; # integra��o num�rica da matriz G
        end
    end
end

h=1/2; # Integra��o anal�tica da matriz H

return g,h
end

function calcula_HeGns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsil,w,k)
    #integra��o n�o singular
    n_pint=length(qsil); # N�mero de pontos de integra��o.
    g=0; # Inicializa o somatorio de g
    h=0; # Inicializa o somatorio de h

    eta = qsil;
    rho = w;
    Tast = complex(0,0);
    qast = complex(0,0);

    for l=1:n_pint # La�o sobre os pontos de integra��o
        for m=1:n_pint # La�o sobre os pontos de integra��o

            qsi=(1-eta[l])*qsil[m];

            N1,N2,N3=calc_fformatri(qsi,eta[l]); #  fun��es de forma
            x=N1*x1+N2*x2+N3*x3; # coordenada x do ponto de integra��o
            y=N1*y1+N2*y2+N3*y3; # coordenada y do ponto de integra��o
            z=N1*z1+N2*z2+N3*z3; # coordenada z do ponto de integra��o

            dNdqsi = [1; 0; -1]; # Derivadas das fun��es de forma
            dNdeta = [0; 1; -1];

            dxdqsi = x1*dNdqsi[1]+x2*dNdqsi[2]+x3*dNdqsi[3];
            dydqsi = y1*dNdqsi[1]+y2*dNdqsi[2]+y3*dNdqsi[3];
            dzdqsi = z1*dNdqsi[1]+z2*dNdqsi[2]+z3*dNdqsi[3];

            dxdeta = x1*dNdeta[1]+x2*dNdeta[2]+x3*dNdeta[3];
            dydeta = y1*dNdeta[1]+y2*dNdeta[2]+y3*dNdeta[3];
            dzdeta = z1*dNdeta[1]+z2*dNdeta[2]+z3*dNdeta[3];

            g1 = dydqsi*dzdeta - dzdqsi*dydeta;
            g2 = dzdqsi*dxdeta - dxdqsi*dzdeta;
            g3 = dxdqsi*dydeta - dydqsi*dxdeta;
            J = sqrt(g1^2.0 + g2^2.0 + g3^2.0);

            Tast,qast=calc_solfund_POT(x,y,z,xd,yd,zd,n,k); # Solu��es fundamentais

            h=h+qast*(1-eta[l])*rho[l]*w[m]*J; # Integral da matriz H
            g=g+Tast*(1-eta[l])*rho[l]*w[m]*J; # Integral da matriz G

        end
    end
    return g,h
end

function calcula_HeGs_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsi,w,k)

# Integra��o singular das matrizes H e G. No caso da matriz G, o elemento
# triangular � decomposto em tr�s quadril�teros degenerados na forma
# de tri�ngulos. Os dois primeiros n�s destes quadril�teros s�o
# coincidentes (formam um s� v�rtice do tri�ngulo) e correspondem ao ponto
# onde existe a singularidade, ou seja, ao ponto fonte que se localiza no
# centr�ide do elemento. Isto faz com que haja uma concentra��o de pontos
# de integra��o junto � singularidade, al�m do jacobiano ser igual a zero
# na singularidade. No caso da matriz H, a integra��o � anal�tica e sempre
# ser� igual a -1/2.

npg=length(qsi); # N�mero de pontos de integra��o
g = 0.0; # inicializa��o da matriz G
Tast = complex(0,0);
qast = complex(0,0);
for kk=1:3
    x1t=xd; # coordenada x do primeiro n� do quadrilatero desgenerado
    y1t=yd; # coordenada y do primeiro n� do quadrilatero desgenerado
    z1t=zd; # coordenada z do primeiro n� do quadrilatero desgenerado
    x2t=xd; # coordenada x do segundo n� do quadrilatero desgenerado
    y2t=yd; # coordenada y do segundo n� do quadrilatero desgenerado
    z2t=zd; # coordenada z do segundo n� do quadrilatero desgenerado

    if(kk==1) # Terceiro e quarto n�s do primeiro quadril�tero desgenerado
        x3t=x1; # coordenada x do terceiro n� do quadrilatero desgenerado
        y3t=y1; # coordenada y do terceiro n� do quadrilatero desgenerado
        z3t=z1; # coordenada z do terceiro n� do quadrilatero desgenerado
        x4t=x2; # coordenada x do quarto n� do quadrilatero desgenerado
        y4t=y2; # coordenada y do quarto n� do quadrilatero desgenerado
        z4t=z2; # coordenada z do quarto n� do quadrilatero desgenerado

    elseif(kk==2) # Terceiro e quarto n�s do segundo quadril�tero
                                  # desgenerado
        x3t=x2; # coordenada x do terceiro n� do quadrilatero desgenerado
        y3t=y2; # coordenada y do terceiro n� do quadrilatero desgenerado
        z3t=z2; # coordenada z do terceiro n� do quadrilatero desgenerado
        x4t=x3; # coordenada x do quarto n� do quadrilatero desgenerado
        y4t=y3; # coordenada y do quarto n� do quadrilatero desgenerado
        z4t=z3; # coordenada z do quarto n� do quadrilatero desgenerado

    elseif(kk==3) # Terceiro e quarto n�s do terceiro quadril�tero
                                  # desgenerado
        x3t=x3; # coordenada x do terceiro n� do quadrilatero desgenerado
        y3t=y3; # coordenada y do terceiro n� do quadrilatero desgenerado
        z3t=z3; # coordenada z do terceiro n� do quadrilatero desgenerado
        x4t=x1; # coordenada x do quarto n� do quadrilatero desgenerado
        y4t=y1; # coordenada y do quarto n� do quadrilatero desgenerado
        z4t=z1; # coordenada z do quarto n� do quadrilatero desgenerado

    end

    for ii = 1: npg; # la�o sobre a primeira vari�vel de integra��o
        for jj = 1: npg # la�o sobre a segunda vari�vel de integra��o
            N1,N2,N3,N4 = calc_fforma(qsi[ii],qsi[jj]); # Fun��es de
                                                       #  forma
            xc = x1t*N1+x2t*N2+x3t*N3+x4t*N4; # coordenada x do
            # ponto campo
            yc = y1t*N1+y2t*N2+y3t*N3+y4t*N4; # coordenada y do
            # ponto campo
            zc = z1t*N1+z2t*N2+z3t*N3+z4t*N4; # coordenada z do
            # ponto campo
            J = calc_jacobiano_quad(x1t,y1t,z1t,x2t,y2t,z2t,x3t,y3t,z3t,x4t,y4t,z4t,qsi[ii],qsi[jj]); # jacobiano(varia ao longo
                         #  do elemento desgenerado)
            Tast,qast = calc_solfund_POT(xd,yd,zd, xc, yc, zc, [0 0 0], k); # Sol.
            # fudamental de temperatura
            g = g + w[jj]*w[ii]*J*Tast; # integra��o num�rica da matriz G
        end
    end
end

h=-1/2; # Integra��o anal�tica da matriz H

return g,h
end

function calcula_Gs(X1,X2,X3,Xd,W,FR,CW,N1,N2,N3,N4,dN1dqsi,dN2dqsi,dN3dqsi,dN4dqsi,dN1deta,dN2deta,dN3deta,dN4deta)
	#Essa eh a funcao que calcula a integral singular da matriz G
	g = 0.0; # inicialização da matriz G

	for kk=1:3
	    x1t=Xd[1]; # coordenada x do primeiro nó do quadrilatero desgenerado
	    y1t=Xd[2]; # coordenada y do primeiro nó do quadrilatero desgenerado
	    z1t=Xd[3]; # coordenada z do primeiro nó do quadrilatero desgenerado
	    x2t=Xd[1]; # coordenada x do segundo nó do quadrilatero desgenerado
	    y2t=Xd[2]; # coordenada y do segundo nó do quadrilatero desgenerado
	    z2t=Xd[3]; # coordenada z do segundo nó do quadrilatero desgenerado

	    if(kk==1) # Terceiro e quarto nós do primeiro quadrilátero desgenerado
	        x3t=X1[1]; # coordenada x do terceiro nó do quadrilatero desgenerado
	        y3t=X1[2]; # coordenada y do terceiro nó do quadrilatero desgenerado
	        z3t=X1[3]; # coordenada z do terceiro nó do quadrilatero desgenerado
	        x4t=X2[1]; # coordenada x do quarto nó do quadrilatero desgenerado
	        y4t=X2[2]; # coordenada y do quarto nó do quadrilatero desgenerado
	        z4t=X2[3]; # coordenada z do quarto nó do quadrilatero desgenerado

	    elseif(kk==2) # Terceiro e quarto nós do segundo quadrilátero
	                                  # desgenerado
	        x3t=X2[1]; # coordenada x do terceiro nó do quadrilatero desgenerado
	        y3t=X2[2]; # coordenada y do terceiro nó do quadrilatero desgenerado
	        z3t=X2[3]; # coordenada z do terceiro nó do quadrilatero desgenerado
	        x4t=X3[1]; # coordenada x do quarto nó do quadrilatero desgenerado
	        y4t=X3[2]; # coordenada y do quarto nó do quadrilatero desgenerado
	        z4t=X3[3]; # coordenada z do quarto nó do quadrilatero desgenerado

	    elseif(kk==3) # Terceiro e quarto nós do terceiro quadrilátero
	                                  # desgenerado
	        x3t=X3[1]; # coordenada x do terceiro nó do quadrilatero desgenerado
	        y3t=X3[2]; # coordenada y do terceiro nó do quadrilatero desgenerado
	        z3t=X3[3]; # coordenada z do terceiro nó do quadrilatero desgenerado
	        x4t=X1[1]; # coordenada x do quarto nó do quadrilatero desgenerado
	        y4t=X1[2]; # coordenada y do quarto nó do quadrilatero desgenerado
	        z4t=X1[3]; # coordenada z do quarto nó do quadrilatero desgenerado

	    end

	    X1t = [x1t, y1t, z1t];
	    X2t = [x2t, y2t, z2t];
	    X3t = [x3t, y3t, z3t];
	    X4t = [x4t, y4t, z4t];

	    x = N1*X1t[1] + N2*X2t[1] + N3*X3t[1] + N4*X4t[1];
	    y = N1*X1t[2] + N2*X2t[2] + N3*X3t[2] + N4*X4t[2];
	    z = N1*X1t[3] + N2*X2t[3] + N3*X3t[3] + N4*X4t[3];


	    # Jacobiano
			dxdqsi = X1t[1]*dN1dqsi + X2t[1]*dN2dqsi + X3t[1]*dN3dqsi + X4t[1]*dN4dqsi;
	    dydqsi = X1t[2]*dN1dqsi + X2t[2]*dN2dqsi + X3t[2]*dN3dqsi + X4t[2]*dN4dqsi;
	    dzdqsi = X1t[3]*dN1dqsi + X2t[3]*dN2dqsi + X3t[3]*dN3dqsi + X4t[3]*dN4dqsi;

	    dxdeta = X1t[1]*dN1deta + X2t[1]*dN2deta + X3t[1]*dN3deta + X4t[1]*dN4deta;
	    dydeta = X1t[2]*dN1deta + X2t[2]*dN2deta + X3t[2]*dN3deta + X4t[2]*dN4deta;
	    dzdeta = X1t[3]*dN1deta + X2t[3]*dN2deta + X3t[3]*dN3deta + X4t[3]*dN4deta;

	    g1 = dydqsi.*dzdeta - dzdqsi.*dydeta;
	    g2 = dzdqsi.*dxdeta - dxdqsi.*dzdeta;
	    g3 = dxdqsi.*dydeta - dydqsi.*dxdeta;
	    J = real(sqrt.(g1.^2.0 + g2.^2.0 + g3.^2.0));

	    # Solução fundamental

	    r1 = x - Xd[1];
	    r2 = y - Xd[2];
	    r3 = z - Xd[3];
	    R = real(sqrt.(r1.^2 + r2.^2 + r3.^2));
#			ZW=complex(0.,-FR*R/CW);
#			past=1; #exp(ZW)/complex(4*pi*R,0.); #Solucao fundamental para a pressao acustica
			Tast = 1.0./(4.0*FR*pi*R); #Solucao fundamental para a temperatura
			g = g + sum(sum(Tast.*W.*J));

	end
return g
end

function calcula_SGeSHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k)
# Evaluates the hypersingular integral kernel
n_pint=length(qsi); # Number of integration points
gx=complex(0,0); # Starts the sum of g in the x direction
gy=complex(0,0); # Starts the sum of g in the y direction
gz=complex(0,0); # Starts the sum of g in the z direction
hx=complex(0,0); # Starts the sum of h in the x direction
hy=complex(0,0); # Starts the sum of h in the y direction
hz=complex(0,0); # Starts the sum of h in the z direction
N = zeros(3)
for l=1:n_pint # Loop over the integration points
    for m=1:n_pint # Loop over the integration points
	xi = (1 - qsi[l])*qsi[m]
        N =calc_fformatri(xi,qsi[l]); #  fun��es de forma
	#dN = calc_dfforma(xi,qsi[m])
        x=N[1]*x1+N[2]*x2+N[3]*x3; # coordenada x do ponto de integra��o
        y=N[1]*y1+N[2]*y2+N[3]*y3 # coordenada y do ponto de integra��o
        z=N[1]*z1+N[2]*z2+N[3]*z3; # coordenada z do ponto de integra��o
	J = calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,xi,qsi[m]);# jacobiano
        #Tast,qast=calc_solfund(x,y,z,xd,yd,zd,n,FR,CW); # Solu��es
        #  fundamentais
	phiastx,phiasty,phiastz,qastx,qasty,qastz=calc_dsolfund(x,y,z,xd,yd,zd,n,k); # Solu��es
        #  fundamentais
	
	#Tast = 1 #To test the numerical integration, the value of the fundamental solutions is set to 1. In this special case, the value of the sum of the rows of matrices H or G should equal the surface area of the model.
	#qast = 1        
	gx=gx+phiastx*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz G
	gy=gy+phiasty*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz G
	gz=gz+phiastz*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz G
        hx=hx+qastx*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz H
        hy=hy+qasty*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz H
        hz=hz+qastz*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz H
    end
end
return gx,gy,gz,hx,hy,hz
end

function calcula_SGeSHns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k)
#integra��o n�o singular
n_pint=length(qsi); # N�mero de pontos de integra��o.
g=complex(0,0); # Inicializa o somatorio de g
h=complex(0,0); # Inicializa o somatorio de h
N = zeros(3)
for l=1:n_pint # La�o sobre os pontos de integra��o
    for m=1:n_pint # La�o sobre os pontos de integra��o
	xi = (1 - qsi[l])*qsi[m]
        N =calc_fformatri(xi,qsi[l]); #  fun��es de forma
	#dN = calc_dfforma(xi,qsi[m])
        x=N[1]*x1+N[2]*x2+N[3]*x3; # coordenada x do ponto de integra��o
        y=N[1]*y1+N[2]*y2+N[3]*y3 # coordenada y do ponto de integra��o
        z=N[1]*z1+N[2]*z2+N[3]*z3; # coordenada z do ponto de integra��o
	J = calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,xi,qsi[m]);# jacobiano
        #Tast,qast=calc_solfund(x,y,z,xd,yd,zd,n,FR,CW); # Solu��es
        #  fundamentais
	phiast,qast=calc_dsolfund_POT(x,y,z,xd,yd,zd,n,k); # Solu��es
        #  fundamentais
	
	#Tast = 1 #To test the numerical integration, the value of the fundamental solutions is set to 1. In this special case, the value of the sum of the rows of matrices H or G should equal the surface area of the model.
	#qast = 1        
	g=g+phiast*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz G
        h=h+qast*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz H
    end
end
return g,h
end

function calc_solfund(x,y,z,xd,yd,zd,n,k)
# Evaluates the fundamental solutions of the Helmholtz equation.

# Determine the distance between source and field points
rx=x-xd;
ry=y-yd;
rz=z-zd;
r =sqrt(rx^2+ry^2+rz^2);

ZW=complex(0.,k*r);
U=exp(ZW)/complex(4*pi*r,0.); # Fundamental solution for the velocity potential
drdn=(rx*n[1] + ry*n[2] + rz*n[3])/r;	# Unit normal derivative of the distance
Q=(ZW-complex(1.,0.))*U*complex(drdn/r,0.); # Fundamental solution for the flux
return U,Q
end

function calc_solfund_POT(x,y,z,xd,yd,zd,n,k)
rx=x-xd;
ry=y-yd;
rz=z-zd;

r =sqrt(rx^2+ry^2+rz^2);
Tast = 1.0/(4.0*pi*r);
qast = (rx*n[1] + ry*n[2] + rz*n[3])/(4.0*pi*r^3.0);

return Tast, qast
end

function calc_q_pint(PONTOS_int,NOS_GEO,ELEM,T,q,k,qsi,w,inc)
# Evaluates the flux of the velocity potential at domain points
n_pint=length(PONTOS_int[:,1]); # Number of domain points
n_elem=length(T); # Number of elements
Gx_int = complex(zeros(n_pint,n_elem))
Gy_int = complex(zeros(n_pint,n_elem))
Gz_int = complex(zeros(n_pint,n_elem))
Hx_int = complex(zeros(n_pint,n_elem))
Hy_int = complex(zeros(n_pint,n_elem))
Hz_int = complex(zeros(n_pint,n_elem))
phi_inc = complex(zeros(n_pint,1))
g = complex(zeros(n_pint,1))
for i=1:n_pint # Loop over the integration points
    x_fonte=PONTOS_int[i,2]; # x coordinate of the source point
    y_fonte=PONTOS_int[i,3]; # y coordinate of the source point
    z_fonte=PONTOS_int[i,4]; # z coordinate of the source point
    for j=1:n_elem  # Loop over the elements
        no1=ELEM[j,2]; # First node of the element
        no2=ELEM[j,3]; # Second node of the element
        no3=ELEM[j,4]; # Third node of the element

        x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
        y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1
        z1=NOS_GEO[no1,4]; # Coordenada z do n� geom�trico 1

        x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
        y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2
        z2=NOS_GEO[no2,4]; # Coordenada x do n� geom�trico 2

        x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 3
        y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 3
        z3=NOS_GEO[no3,4]; # Coordenada z do n� geom�trico 3

        n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio
        Gx_int[i,j],Gy_int[i,j],Gz_int[i,j],Hx_int[i,j],Hy_int[i,j],Hz_int[i,j]=calcula_SGeSHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,x_fonte,y_fonte,z_fonte,n,qsi,w,k); # Evaluates the kernel for the hypersingular boundary equation
    end
#    if(inc[1,1]!=0)
#    	#g[i,1]=calc_q(x_fonte,y_fonte,fc,FR,CW,GE);
#			g[i,1]=0;
#		else
#    	g[i,1]=0;
#    end
#		if inc[1,1] != 0
	#Vamos incluir um termo de onda incidente
#				phi_inc[i,1] = calc_inc(x_fonte,y_fonte,z_fonte,FR,CW,inc[1,:]);
#		end
end
# Velocity potential flux at domain points
dphidx=Hx_int*T-Gx_int*q;
dphidy=Hy_int*T-Gy_int*q; 
dphidz=Hz_int*T-Gz_int*q; 
#T_pint = - (H_int*T - G_int*q - phi_inc)
#T_pint=-(H_int*T'-G_int*q'-g'); 
return dphidx,dphidy,dphidz
end

function calc_q_pint_POT(PONTOS_int,NOS_GEO,ELEM,T,q,k,qsi,w,inc)
# Calcula a temperatura nos pontos internos
n_pint=length(PONTOS_int[:,1]); # Numero de pontos internos
n_elem=length(T); # Numero de elementos
G_int = complex(zeros(n_pint,n_elem))
H_int = complex(zeros(n_pint,n_elem))
phi_inc = complex(zeros(n_pint,1))
g = complex(zeros(n_pint,1))
for i=1:n_pint # La�o sobre os pontos internos
    x_fonte=PONTOS_int[i,2]; # Coordenada x do ponto fonte
    y_fonte=PONTOS_int[i,3]; # Coordenada y do ponto fonte
    z_fonte=PONTOS_int[i,4]; # Coordenada z do ponto fonte
    for j=1:n_elem  #La�o sobre os elementos
        no1=ELEM[j,2]; # Primeiro n� geom�trico do elemento
        no2=ELEM[j,3]; # Segundo n� geom�trico do elemento
        no3=ELEM[j,4]; # Terceiro n� geom�trico do elemento

        x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
        y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1
        z1=NOS_GEO[no1,4]; # Coordenada z do n� geom�trico 1

        x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
        y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2
        z2=NOS_GEO[no2,4]; # Coordenada x do n� geom�trico 2

        x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 3
        y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 3
        z3=NOS_GEO[no3,4]; # Coordenada z do n� geom�trico 3

        n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio
        G_int[i,j],H_int[i,j]=calcula_SGeSHns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,x_fonte,y_fonte,z_fonte,n,qsi,w,k); # Chama a functio para calculo de H e G
        # quando o ponto fonte nao pertence ao elemento
    end
#    if(inc[1,1]!=0)
#    	#g[i,1]=calc_q(x_fonte,y_fonte,fc,FR,CW,GE);
#			g[i,1]=0;
#		else
#    	g[i,1]=0;
#    end
#		if inc[1,1] != 0
	#Vamos incluir um termo de onda incidente
#				phi_inc[i,1] = calc_inc(x_fonte,y_fonte,z_fonte,FR,CW,inc[1,:]);
#		end
end
T_pint = (H_int*T - G_int*q - phi_inc)
#T_pint=-(H_int*T'-G_int*q'-g'); # Vetor que contem a temperatura nos
#      pontos internos
return T_pint
end

function calc_T_pint(PONTOS_int,NOS_GEO,ELEM,T,q,k,qsi,w,inc)
# Calcula a temperatura nos pontos internos
n_pint=length(PONTOS_int[:,1]); # Numero de pontos internos
n_elem=length(T); # Numero de elementos
G_int = complex(zeros(n_pint,n_elem))
H_int = complex(zeros(n_pint,n_elem))
phi_inc = complex(zeros(n_pint,1))
g = complex(zeros(n_pint,1))
for i=1:n_pint # La�o sobre os pontos internos
    x_fonte=PONTOS_int[i,2]; # Coordenada x do ponto fonte
    y_fonte=PONTOS_int[i,3]; # Coordenada y do ponto fonte
    z_fonte=PONTOS_int[i,4]; # Coordenada z do ponto fonte
    for j=1:n_elem  #La�o sobre os elementos
        no1=ELEM[j,2]; # Primeiro n� geom�trico do elemento
        no2=ELEM[j,3]; # Segundo n� geom�trico do elemento
        no3=ELEM[j,4]; # Terceiro n� geom�trico do elemento

        x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
        y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1
        z1=NOS_GEO[no1,4]; # Coordenada z do n� geom�trico 1

        x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
        y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2
        z2=NOS_GEO[no2,4]; # Coordenada x do n� geom�trico 2

        x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 3
        y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 3
        z3=NOS_GEO[no3,4]; # Coordenada z do n� geom�trico 3

        n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio
        G_int[i,j],H_int[i,j]=calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,x_fonte,y_fonte,z_fonte,n,qsi,w,k); # Chama a functio para calculo de H e G
        # quando o ponto fonte nao pertence ao elemento
    end
#    if(inc[1,1]!=0)
#    	#g[i,1]=calc_q(x_fonte,y_fonte,fc,FR,CW,GE);
#			g[i,1]=0;
#		else
#    	g[i,1]=0;
#    end
#		if inc[1,1] != 0
	#Vamos incluir um termo de onda incidente
#				phi_inc[i,1] = calc_inc(x_fonte,y_fonte,z_fonte,FR,CW,inc[1,:]);
#		end
end
T_pint = - (H_int*T - G_int*q - phi_inc)
#T_pint=-(H_int*T'-G_int*q'-g'); # Vetor que contem a temperatura nos
#      pontos internos
return T_pint
end

function calc_T_pint_POT(PONTOS_int,NOS_GEO,ELEM,T,q,k,qsi,w,inc)
# Calcula a temperatura nos pontos internos
n_pint=length(PONTOS_int[:,1]); # Numero de pontos internos
n_elem=length(T); # Numero de elementos
G_int = complex(zeros(n_pint,n_elem))
H_int = complex(zeros(n_pint,n_elem))
phi_inc = complex(zeros(n_pint,1))
g = complex(zeros(n_pint,1))
for i=1:n_pint # La�o sobre os pontos internos
    x_fonte=PONTOS_int[i,2]; # Coordenada x do ponto fonte
    y_fonte=PONTOS_int[i,3]; # Coordenada y do ponto fonte
    z_fonte=PONTOS_int[i,4]; # Coordenada z do ponto fonte
    for j=1:n_elem  #La�o sobre os elementos
        no1=ELEM[j,2]; # Primeiro n� geom�trico do elemento
        no2=ELEM[j,3]; # Segundo n� geom�trico do elemento
        no3=ELEM[j,4]; # Terceiro n� geom�trico do elemento

        x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
        y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1
        z1=NOS_GEO[no1,4]; # Coordenada z do n� geom�trico 1

        x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
        y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2
        z2=NOS_GEO[no2,4]; # Coordenada x do n� geom�trico 2

        x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 3
        y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 3
        z3=NOS_GEO[no3,4]; # Coordenada z do n� geom�trico 3

        n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio
        G_int[i,j],H_int[i,j]=calcula_HeGns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,x_fonte,y_fonte,z_fonte,n,qsi,w,k); # Chama a functio para calculo de H e G
        # quando o ponto fonte nao pertence ao elemento
    end
    if(inc[1,1]!=0)
    	#g[i,1]=calc_q(x_fonte,y_fonte,fc,FR,CW,GE);
			g[i,1]=0;
		else
    	g[i,1]=0;
    end
		if inc[1,1] != 0
	#Vamos incluir um termo de onda incidente
				phi_inc[i,1] = calc_inc(x_fonte,y_fonte,z_fonte,FR,CW,inc[1,:]);
		end
end
T_pint =  (H_int*T - G_int*q)
#T_pint=-(H_int*T'-G_int*q'-g'); # Vetor que contem a temperatura nos
#      pontos internos
return T_pint
end

function calc_dsolfund(x,y,z,xd,yd,zd,n,k)
# Evaluates the fundamental solutions for the derivative of the Helmholtz equation.

# Determine the distance between source and field points
rx=x-xd;
ry=y-yd;
rz=z-zd;
r =sqrt(rx^2+ry^2+rz^2);

kr = complex(0,-k*r);
drdn=(rx*n[1] + ry*n[2] + rz*n[3]); # Unit normal derivative of the distance
U=exp(kr)./complex(4*pi*r,0.);
phiast_x = (kr-complex(1.,0.)).*U.*complex(rx./r.^2,0.);
phiast_y = (kr-complex(1.,0.)).*U.*complex(ry./r.^2,0.);
phiast_z = (kr-complex(1.,0.)).*U.*complex(rz./r.^2,0.);
qast_x = (1./(4*pi))*(-n[1]*(1./r.^3+complex(0,1)*k./(r.^2)) + (3./r.^3+3*complex(0,1)*k./(r.^2)-k^2./(r)).*rx.*drdn./r.^2).*exp(-complex(0,1)*k*r);
qast_y = (1./(4*pi))*(-n[2]*(1./r.^3+complex(0,1)*k./(r.^2)) + (3./r.^3+3*complex(0,1)*k./(r.^2)-k^2./(r)).*ry.*drdn./r.^2).*exp(-complex(0,1)*k*r);
qast_z = (1./(4*pi))*(-n[3]*(1./r.^3+complex(0,1)*k./(r.^2)) + (3./r.^3+3*complex(0,1)*k./(r.^2)-k^2./(r)).*rz.*drdn./r.^2).*exp(-complex(0,1)*k*r);
return phiast_x,phiast_y,phiast_z,qast_x,qast_y,qast_z
end

function calc_dsolfund_POT(x,y,z,xd,yd,zd,n,k)
# Evaluates the fundamental solutions for the derivative of the Helmholtz equation.
# Determine the distance between source and field points
rx=x-xd;
ry=y-yd;
rz=z-zd;
r =sqrt(rx^2+ry^2+rz^2);

drdn=(rx*n[1] + ry*n[2] + rz*n[3])/r; # Unit normal derivative of the distance
phiast = (1/(4*pi*r^2))*drdn;
qastx = (1/(4*pi*r^3))*(n[1] - rx*drdn);
qasty = (1/(4*pi*r^3))*(n[2] - ry*drdn);
qastz = (1/(4*pi*r^3))*(n[3] - rz*drdn);
qast = sqrt(qastx^2 +qasty^2 +qastz^2)
return phiast,qast
end

function calc_inc(xd,yd,zd,FR,CW,inc)
#Funcao para calcular o valor da influencia de uma onda incidente no
#elemento estudado
k = complex(0,FR/CW); #Numero de onda
d = [inc[2,1] inc[3,1] inc[4,1]]; #direcao de propagacao da onda incidente /d/ = 1
p = [xd yd zd];
dp = dot(d,p);
A_inc = inc[5,1]; #amplitude da onda
phi_inc = A_inc*exp(k*dp);
return phi_inc
end

function calc_fforma(qsi, eta)

#--------------------------------------------------------------------------
# Dados de entrada:
# qsi, eta - pontos de Gauss onde as fun��es de forma s�o calculadas.
# Dados de sa�da:
# [N] - Fun��es de forma para um elemento quadrilateral linear calculadas
# em (qsi, eta).
#--------------------------------------------------------------------------

N = (1./4.)*[(1. - qsi).*(1. - eta);
               (1. + qsi).*(1. - eta);
               (1. + qsi).*(1. + eta);
               (1. - qsi).*(1. + eta)];
N1=N[1];N2=N[2];N3=N[3];N4=N[4];
return N1,N2,N3,N4
end

function calc_fformatri(qsi,eta)

    N1 = qsi;
    N2 = eta;
    N3 = 1-qsi-eta;

return N1,N2,N3
end

function calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,qsi,eta)
dNdqsi, dNdeta = calc_dfforma(qsi,eta); # Calcula a derivada das fun��es
  # de forma
dxdqsi = x1*dNdqsi[1]+x2*dNdqsi[2]+x3*dNdqsi[3]
dydqsi = y1*dNdqsi[1]+y2*dNdqsi[2]+y3*dNdqsi[3]
dzdqsi = z1*dNdqsi[1]+z2*dNdqsi[2]+z3*dNdqsi[3]

dxdeta = x1*dNdeta[1]+x2*dNdeta[2]+x3*dNdeta[3]
dydeta = y1*dNdeta[1]+y2*dNdeta[2]+y3*dNdeta[3]
dzdeta = z1*dNdeta[1]+z2*dNdeta[2]+z3*dNdeta[3]

g1 = dydqsi*dzdeta - dzdqsi*dydeta;
g2 = dzdqsi*dxdeta - dxdqsi*dzdeta;
g3 = dxdqsi*dydeta - dydqsi*dxdeta;
J = sqrt(g1^2.0 + g2^2.0 + g3^2.0);

return J
end

function calc_jacobiano_quad(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,qsi,eta)
dNdqsi, dNdeta = calc_dfforma_quad(qsi,eta); # Calcula a derivada das fun��es
  # de forma
dxdqsi = x1*dNdqsi[1]+x2*dNdqsi[2]+x3*dNdqsi[3]+x4*dNdqsi[4];
dydqsi = y1*dNdqsi[1]+y2*dNdqsi[2]+y3*dNdqsi[3]+y4*dNdqsi[4];
dzdqsi = z1*dNdqsi[1]+z2*dNdqsi[2]+z3*dNdqsi[3]+z4*dNdqsi[4];

dxdeta = x1*dNdeta[1]+x2*dNdeta[2]+x3*dNdeta[3]+x4*dNdeta[4];
dydeta = y1*dNdeta[1]+y2*dNdeta[2]+y3*dNdeta[3]+y4*dNdeta[4];
dzdeta = z1*dNdeta[1]+z2*dNdeta[2]+z3*dNdeta[3]+z4*dNdeta[4];

g1 = dydqsi*dzdeta - dzdqsi*dydeta;
g2 = dzdqsi*dxdeta - dxdqsi*dzdeta;
g3 = dxdqsi*dydeta - dydqsi*dxdeta;
J = sqrt(g1^2.0 + g2^2.0 + g3^2.0);

return J
end

function calc_dfforma(qsi, eta)

dNdqsi = [1 0 -1]

dNdeta = [0 1 -1]
return dNdqsi, dNdeta
end
function calc_dfforma_quad(qsi, eta)

dNdqsi = (1./4.)*[-(1. - eta);
                  (1. - eta);
                  (1. + eta);
                 -(1. + eta)];

dNdeta = (1./4.)*[-(1. - qsi);
                 -(1. + qsi);
                  (1. + qsi);
                  (1. - qsi)];
return dNdqsi, dNdeta
end

function calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3)
    # Function que calcula o vetor unit�rio normal ao elemento

    v1 = [x3,y3,z3] - [x2,y2,z2]; # vetor formado pela aresta 32 do elemento
    v2 = [x1,y1,z1] - [x2,y2,z2]; # vetor formado pela aresta 12 do elemento
    n = cross(v1, v2); # Produto vetorial entre v1 e v2 (vetor normal ao
    # elemento)
    #println(v1,v2,n)
    if sum(abs.(n)) > 0.000001
	n = n./norm(n); # vetor unit�rio normal ao elemento
    end
    #println(v1,v2,n)
    return n
end

function Gauss_Legendre(x1,x2,n)
eps=3e-14;
m::Int64 = round((n+1)/2);
x = zeros(1,n)
w = zeros(1,n)
pp = 1
xm=0.5*(x2+x1);
xl=0.5*(x2-x1);
for i=1:m
    z=cos(pi*(i-0.25)/(n+0.5));
    while 1==1
        p1=1.0;
        p2=0.0;
        for j=1:n
            p3=p2;
            p2=p1;
            p1=((2*j-1)*z*p2-(j-1)*p3)/j;
        end
        pp=n*(z*p1-p2)/(z*z-1);
        z1=z;
        z=z1-p1/pp;
        if(abs(z-z1)<eps)
            break
        end
    end
    x[i]=xm-xl*z;
    x[n+1-i]=xm+xl*z;
    w[i]=2*xl/((1-z*z)*pp*pp);
    w[n+1-i]=w[i];
end
return x,w
end

function gera_CDC(ELEM,CCFace)

# Gera a matriz de CDC
# CDC = [número do elemento, tipo da CDC, valor da CDC no nó 1,...
#                               valor da CDC no nó 2, valor da CDC no nó 3]
nelem = length(ELEM[:,1]);
CDC = zeros(nelem,3);
for i=ELEM[:,1]
    CDC[i,1] = i;
    CDC[i,2] = CCFace[ELEM[i,5],2];
    CDC[i,3] = CCFace[ELEM[i,5],3];
end
return CDC
end

function mini(A)
min = Inf;
for i = 1:size(A,1)
	if abs(A[i,1])<min
		min = A[i,1];
	end
end
return mini
end

function maxi(A)
maxi = 0.0;
for i = 1:size(A,1)
	if abs(A[i,1])>maxi
		min = A[i,1];
	end
end
return maxi
end

function montamatrizvetor(NOS, NOS_GEO, ELEM, k, CDC)
    nelem=size(ELEM,1); # N�mero de elementos (n�mero de linhas da
    #  matriz ELEM)

    nnos=nelem; # N�mero de n�s
    G=zeros(nnos,nnos); # Inicializa��o da matriz G
    H=zeros(nnos,nnos); # Inicializa��o da matriz H

    npg=12; # N�mero de pontos de Gauss
    qsi_tri,w_tri=Gauss_Legendre(0,1,npg); # Pontos e pesos de Gauss do
    # elemento triangular
    qsi,w=Gauss_Legendre(-1,1,npg); # Pontos e pesos de Gauss do elemento
    # quadrilateral
    A = complex(zeros(nnos,nnos));
    b = complex(zeros(nnos, 1));

    for j=1:nelem # La�o sobre os elementos
	tipoCDC = CDC[j,2]; #Tipo da condicao de contorno CDC[pc,2] = 0 condicao de pressao, =1 condicao de fluxo
	#valorCDC = complex(CDC[pc,3],CDC[pc,4])	#Valor da condicao no elementos
	valorCDC = CDC[j,3];
	nos = ELEM[j,2:4];		#Nos contidos no elementos

        no1=ELEM[j,2]; # Primeiro n� geom�trico do elemento
        no2=ELEM[j,3]; # Segundo n� geom�trico do elemento
        no3=ELEM[j,4]; # Terceiro n� geom�trico do elemento

        x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
        y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1
        z1=NOS_GEO[no1,4]; # Coordenada z do n� geom�trico 1

        x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
        y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2
        z2=NOS_GEO[no2,4]; # Coordenada x do n� geom�trico 2

        x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 3
        y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 3
        z3=NOS_GEO[no3,4]; # Coordenada z do n� geom�trico 3

        n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio
        # normal ao elemento
        for i=1:nnos # La�o sobre os pontos fontes
            xd=NOS[i,2]; # Coordenada x do ponto fonte
            yd=NOS[i,3]; # Coordenada y do ponto fonte
            zd=NOS[i,4]; # Coordenada y do ponto fonte

            if i==j # O ponto fonte pertence ao elemento
                G[i,j],H[i,j]=calcula_HeGs_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsi,w,k); # Integra��o singular
            else # O ponto fonte n�o pertence ao elemento
                G[i,j],H[i,j]=calcula_HeGns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi_tri,w_tri,k); # Integra��o regular
            end
            #Aplica as condicoes de contorno
            if tipoCDC == 0 # A pressao eh conhecida
        	A[i,j] = -G[i,j]; 	#Os valores de G vao para a matriz A
                b[i,1] = b[i,1] - H[i,j]*valorCDC; #Os valores de H vao para o vetor b
            else
                A[i,j] = +H[i,j]; 	#Os valores de H vao para a matriz A
                b[i,1] = b[i,1] + G[i,j]*valorCDC; #Os valores de G vao para o vetor b
            end
        end
    end
    return G,H
end

function monta_matrizvetor(NOS,NOS_GEO, ELEM, FR, CW, qsi, w, qsi_quad,w_quad, inc)
nelem=size(ELEM,1); #Numero de elementos
nnos = nelem; #No caso de elementos constantes, o numero de nos eh o mesmo que o de elementos
A = complex.(zeros(nnos,nnos));
b = complex.(zeros(nnos, 1));
#Inicia os pontos e pesos de Gauss para os elementos triangulares lineares
W = w'*w;
eta = qsi;
ETA = ones(size(qsi))'*qsi;
#Calcula as funcoes de forma triangulares
N1_tri = ((1-eta)'*qsi)';
N2_tri = ones(length(w),1)*eta;
N3_tri = 1 - N1_tri - N2_tri;
#Para a integracao singular, vamos desgenerar em um elemento quadrilateral
W_quad = w_quad'*w_quad;
#Calcula as funcoes de forma quadrilaterais
N1 = (1./4.)*(1. - qsi_quad')*(1. - qsi_quad);
N2=(1./4.)*(1. + qsi_quad')*(1. - qsi_quad);
N3=(1./4.)*(1. + qsi_quad')*(1. + qsi_quad);
N4=(1./4.)*(1. - qsi_quad')*(1. + qsi_quad);

dN1dqsi = (1./4.)*ones(npg,1)*(-(1. - qsi_quad));
dN2dqsi = (1./4.)*ones(npg,1)*(1. - qsi_quad);
dN3dqsi = (1./4.)*ones(npg,1)*(1. + qsi_quad);
dN4dqsi = (1./4.)*ones(npg,1)*(-(1. + qsi_quad));

dN1deta = (1./4.)*ones(npg,1)*(-(1. - qsi_quad));
dN2deta = (1./4.)*ones(npg,1)*(-(1. + qsi_quad));
dN3deta = (1./4.)*ones(npg,1)*(1. + qsi_quad);
dN4deta = (1./4.)*ones(npg,1)*(1. - qsi_quad);

# Inicio da integracao

for pc = 1:nelem #Laco sobre os elementos
#println("pc= ",pc)
	tipoCDC = CDC[pc,2] #Tipo da condicao de contorno CDC[pc,2] = 0 condicao de pressao, =1 condicao de fluxo
	#valorCDC = complex(CDC[pc,3],CDC[pc,4])	#Valor da condicao no elementos
	valorCDC = CDC[pc,3];
	nos = ELEM[pc,2:4]		#Nos contidos no elementos

	X1 = NOS_GEO[nos[1],2:4] 	#Coordenada do primeiro no
	X2 = NOS_GEO[nos[2],2:4] 	#Coordenada do primeiro no
	X3 = NOS_GEO[nos[3],2:4] 	#Coordenada do primeiro no
	v1 = X3 - X2; 	#vetor formado pela aresta 3-2 do elemento
	v2 = X1 - X2; #vetor formado pela aresta 1-2 do elemento
	n = [v1[2]*v2[3] - v2[2]*v1[3]; v2[1]*v1[3]-v1[1]*v2[3]; v1[1]*v2[2] - v2[1]*v1[2]];
	n = n./norm(n);
	J =     J=real(sqrt((-(X1[2]*X2[1])+X1[1]*X2[2]+X1[2]*X3[1]-X2[2]*X3[1] - X1[1]*X3[2]+X2[1]*X3[2])^2+(X1[3]*X2[1]-X1[1]*X2[3]-X1[3]*X3[1]+ X2[3]*X3[1]+X1[1]*X3[3]-X2[1]*X3[3])^2+(-(X1[3]*X2[2])+ X1[2]*X2[3]+X1[3]*X3[2]-X2[3]*X3[2]-X1[2]*X3[3]+X2[2]*X3[3])^2));
	x = N1_tri*X1[1] + N2_tri*X2[1] + N3_tri*X3[1]; #Coordenadas x para o elemento
	y = N1_tri*X1[2] + N2_tri*X2[2] + N3_tri*X3[2]; #Coordenadas x para o elemento
	z = N1_tri*X1[3] + N2_tri*X2[3] + N3_tri*X3[3]; #Coordenadas x para o elemento

	g = Array{Any}
	h = Array{Any}
	for pf = 1: nnos 	#Laco sobre os nos
#println("pf= ",pf)
		Xd = NOS[pf,2:4]; #Coordenadas do ponto fonte
		if pf == pc # Para elementos constantes, se o no corresponder ao elemento, a integracao eh singular
			g = calcula_Gs(X1,X2,X3,Xd,W_quad,FR,CW,N1,N2,N3,N4,dN1dqsi,dN2dqsi,dN3dqsi,dN4dqsi,dN1deta,dN2deta,dN3deta,dN4deta);
			h = -1/2;
		#println("singular")
		else
		#println("não singular")
			r1 = x - Xd[1];
			r2 = y - Xd[2];
			r3 = z - Xd[3];
			R = real(sqrt.(r1.^2 + r2.^2 + r3.^2));
			ZW=complex.(0.,-FR*R/CW);
#			past=1; #exp(ZW)./complex(4*pi*R,0.);
#			drdn=(r1*n[1] + r2*n[2] + r3*n[3])./R;
#			qast=1; #(ZW-complex(1.,0.))*past*complex(drdn./R,0.);	#Solucao fundamental para o fluxo da pressao acustica
			Tast = 1.0./(4.0*FR*pi*R); #Solucao fundamental para a temperatura
			qast = (r1*n[1] + r2*n[2] + r3*n[3])./(4.0*pi*R.^3); #Solucao fundamental para o fluxo da temperatura
			g = sum(sum(Tast.*(1-ETA).*W.*J));
			h = sum(sum(qast.*(1-ETA).*W.*J));
		end
	if tipoCDC == 0 # A pressao eh conhecida
		A[pf,pc] = -g; 	#Os valores de G vao para a matriz A
		b[pf,1] = b[pf,1] - h*valorCDC; #Os valores de H vao para o vetor b
#println("CDC = 0. A[",pf,",",pc,"]= ",A[pf,pc])
#println("CDC = 0. b[",pf,",1]= ",b[pf,1])
	else
			A[pf,pc] = +h; 	#Os valores de H vao para a matriz A
			b[pf,1] = b[pf,1] + g*valorCDC; #Os valores de G vao para o vetor b
	#		println("CDC = 1. A[",pf,",",pc,"]= ",A[pf,pc])
		#	println("CDC = 1. b[",pf,",1]= ",b[pf,1])
	end
end
end
return A,b
end

function monta_Teq(CDC,x)
# Separa fluxo e temperatura

# ncdc = n�mero de linhas da matriz CDC
# T = vetor que cont�m as temperaturas nos n�s
# q = vetor que cont�m o fluxo nos n�s

ncdc = length(CDC[:,1]);
nnos = ncdc
T = complex(zeros(nnos,1))
q = complex(zeros(nnos,1))
for i=1:ncdc # La�o sobre as condi��es de contorno
    tipoCDC=CDC[i,2]; # Tipo da condi��o de contorno
    #valorCDC=complex(CDC[i,3],CDC[i,4]); # Valor da condi��o de contorno
		valorCDC = CDC[i,3];
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

function cal_GeH(NOS,NOS_GEO,ELEM,k,qsi,w,inc)
# Evaluates the G and H matrices for the linear system H phi = G q, where phi is a vector containing the values of the velocity potential and q is a vector containing the values of the flux at the boundary of the problem.

nelem=size(ELEM,1); # N�mero de elementos (n�mero de linhas da
#  matriz ELEM)
qsitelles,Jtelles = telles(qsi,0); # Evaluates the Telles' points and Jacobian
npg=12; # Number of integration points
qsiquad,wquad = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
nnos=nelem; # N�mero de n�s
G=complex(zeros(nnos,nnos)); # Inicializa��o da matriz G
H=complex(zeros(nnos,nnos)); # Inicializa��o da matriz H
phi_inc = complex(zeros(nelem,1));
for i=1:nnos # La�o sobre os pontos fontes
		xd=NOS[i,2]; # Coordenada x do ponto fonte
		yd=NOS[i,3]; # Coordenada y do ponto fonte
		zd=NOS[i,4]; # Coordenada y do ponto fonte

		for j=1:nelem # La�o sobre os elementos
		    no1=ELEM[j,2]; # Primeiro n� geom�trico do elemento
		    no2=ELEM[j,3]; # Segundo n� geom�trico do elemento
		    no3=ELEM[j,4]; # Terceiro n� geom�trico do elemento

		    x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
		    y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1
		    z1=NOS_GEO[no1,4]; # Coordenada z do n� geom�trico 1

		    x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
		    y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2
		    z2=NOS_GEO[no2,4]; # Coordenada x do n� geom�trico 2

		    x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 3
		    y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 3
		    z3=NOS_GEO[no3,4]; # Coordenada z do n� geom�trico 3

		    n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio normal ao elemento
		        if i==j # O ponto fonte pertence ao elemento
		           #G[i,j],H[i,j]=calcula_HeGs_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsiquad,wquad,k); # Integra��o singular
			   G[i,j],H[i,j]= calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsiquad,wquad,k)# Integra��o singular
			   #G[i,j],H[i,j]=calcula_GeHns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k); # Integra��o singular
			    #Gtelles,Htelles=calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsitelles,w.*Jtelles,k); # Integra��o singular
#G[i,j]=1
#H[i,j]= -0.5
#			erro = (G[i,j] - Gtelles)
#			println("diferença entre g e gtelles= $erro")
		        else # O ponto fonte n�o pertence ao elemento
		           # G[i,j],H[i,j]=calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k); # Integra��o
		            #  regular
			   G[i,j],H[i,j]=calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k); # Integra��o singular
		        end
		    end
				if inc[1,1] != 0
			#Vamos incluir um termo de onda incidente
					phi_inc[i,1] = calc_inc(xd,yd,zd,k,k,inc[1,:]);
				end
end

return G,H,phi_inc
end

function cal_GeH_POT(NOS,NOS_GEO,ELEM,k,qsi,w,inc)
# Evaluates the G and H matrices for the linear system H phi = G q, where phi is a vector containing the values of the velocity potential and q is a vector containing the values of the flux at the boundary of the problem.

nelem=size(ELEM,1); # N�mero de elementos (n�mero de linhas da
#  matriz ELEM)
qsitelles,Jtelles = telles(qsi,0); # Evaluates the Telles' points and Jacobian
npg=4; # Number of integration points
qsiquad,wquad = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
nnos=nelem; # N�mero de n�s
G=complex(zeros(nnos,nnos)); # Inicializa��o da matriz G
H=complex(zeros(nnos,nnos)); # Inicializa��o da matriz H
phi_inc = complex(zeros(nelem,1));
for i=1:nnos # La�o sobre os pontos fontes
		xd=NOS[i,2]; # Coordenada x do ponto fonte
		yd=NOS[i,3]; # Coordenada y do ponto fonte
		zd=NOS[i,4]; # Coordenada y do ponto fonte

		for j=1:nelem # La�o sobre os elementos
		    no1=ELEM[j,2]; # Primeiro n� geom�trico do elemento
		    no2=ELEM[j,3]; # Segundo n� geom�trico do elemento
		    no3=ELEM[j,4]; # Terceiro n� geom�trico do elemento

		    x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
		    y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1
		    z1=NOS_GEO[no1,4]; # Coordenada z do n� geom�trico 1

		    x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
		    y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2
		    z2=NOS_GEO[no2,4]; # Coordenada x do n� geom�trico 2

		    x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 3
		    y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 3
		    z3=NOS_GEO[no3,4]; # Coordenada z do n� geom�trico 3

		    n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio normal ao elemento
		        if i==j # O ponto fonte pertence ao elemento
		           #G[i,j],H[i,j]=calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsiquad,wquad,FR,CW); # Integra��o singular
			   G[i,j],H[i,j]= calcula_HeGs_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsiquad,wquad,k)# Integra��o singular
			    #G[i,j],H[i,j]=calcula_GeHns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k); # Integra��o singular
			    #Gtelles,Htelles=calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsitelles,w.*Jtelles,k); # Integra��o singular
#G[i,j]=1
#H[i,j]= -0.5
#			erro = (G[i,j] - Gtelles)
#			println("diferença entre g e gtelles= $erro")
		        else # O ponto fonte n�o pertence ao elemento
		            G[i,j],H[i,j]=calcula_HeGns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k); # Integra��o
		            #  regular
		        end
		    end
				if inc[1,1] != 0
			#Vamos incluir um termo de onda incidente
					phi_inc[i,1] = calc_inc(xd,yd,zd,FR,CW,inc[1,:]);
				end
end

return G,H,phi_inc
end

