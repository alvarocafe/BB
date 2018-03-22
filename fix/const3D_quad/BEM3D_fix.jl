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
return A,b
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
            J = calc_jacobiano(x1t,y1t,z1t,x2t,y2t,z2t,x3t,y3t,z3t,x4t,y4t,z4t,qsi[ii],qsi[jj]); # jacobiano(varia ao longo
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
function calc_solfund_POT(x,y,z,xd,yd,zd,n,k)
rx=x-xd;
ry=y-yd;
rz=z-zd;

r =sqrt(rx^2+ry^2+rz^2);
Tast = 1.0/(4.0*k*pi*r);
qast = (rx*n[1] + ry*n[2] + rz*n[3])/(4.0*pi*r^3.0);

return Tast, qast
end


function calc_fformatri(qsi,eta)

    N1 = qsi;
    N2 = eta;
    N3 = 1-qsi-eta;

return N1,N2,N3
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

function gera_CDC(ELEM,CCFace)

# Gera a matriz de CDC
# CDC = [número do elemento, tipo da CDC, valor da CDC no nó 1,...
#                               valor da CDC no nó 2, valor da CDC no nó 3]
nelem = length(ELEM[:,1]);
CDC = zeros(nelem,3);
for i=ELEM[:,1]
    CDC[i,1] = i;
    CDC[i,2] = CCFace[ELEM[i,6],2];
    CDC[i,3] = CCFace[ELEM[i,6],3];
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

function mostra_geo(NOS_GEO,ELEM)
# geometria da estrutura

nelem=size(ELEM,1);
xmin = mini(NOS_GEO[:,2]);
xmax = maxi(NOS_GEO[:,2]);
ymin = mini(NOS_GEO[:,3]);
ymax = maxi(NOS_GEO[:,3]);
zmin = mini(NOS_GEO[:,4]);
zmax = maxi(NOS_GEO[:,4]);
NOS=zeros(nelem,4);

for i = 1:nelem
    no1=ELEM[i,2];
    no2=ELEM[i,3];
    no3=ELEM[i,4];
    no4=ELEM[i,5];

    xno1=NOS_GEO[no1,2];
    yno1=NOS_GEO[no1,3];
    zno1=NOS_GEO[no1,4];

    xno2=NOS_GEO[no2,2];
    yno2=NOS_GEO[no2,3];
    zno2=NOS_GEO[no2,4];

    xno3=NOS_GEO[no3,2];
    yno3=NOS_GEO[no3,3];
    zno3=NOS_GEO[no3,4];

    xno4=NOS_GEO[no4,2];
    yno4=NOS_GEO[no4,3];
    zno4=NOS_GEO[no4,4];

    xnos=[xno1,xno2,xno3,xno4];
    ynos=[yno1,yno2,yno3,yno4];
    znos=[zno1,zno2,zno3,zno4];
    #xlabel ('{\it x}')
    #ylabel ('{\it y}')
    #zlabel ('{\it z}')
    #fill3(xnos,ynos,znos,'w');	#Parte do programa que mostra geometria
    #axis equal
    #axis([xmin xmax ymin ymax zmin zmax])
    #hold on;
    xno=sum(xnos)/4;
    yno=sum(ynos)/4;
    zno=sum(znos)/4;
    #text(xno,yno,zno,num2str(i))
    NOS[i,:]=[i,xno,yno,zno];
end;
#view(3)
#hold off;

return NOS
end

function mostra_geoTRI(NOS_GEO,ELEM)
# geometria da estrutura

nelem=size(ELEM,1);
# xmin = min(NOS_GEO(:,2));
# xmax = max(NOS_GEO(:,2));
# ymin = min(NOS_GEO(:,3));
# ymax = max(NOS_GEO(:,3));
# zmin = min(NOS_GEO(:,4));
# zmax = max(NOS_GEO(:,4));
NOS=zeros(nelem,4);

for i = 1:nelem
		no1=ELEM[i,2];
		no2=ELEM[i,3];
		no3=ELEM[i,4];

		xno1=NOS_GEO[no1,2];
		yno1=NOS_GEO[no1,3];
		zno1=NOS_GEO[no1,4];

		xno2=NOS_GEO[no2,2];
		yno2=NOS_GEO[no2,3];
		zno2=NOS_GEO[no2,4];

		xno3=NOS_GEO[no3,2];
		yno3=NOS_GEO[no3,3];
		zno3=NOS_GEO[no3,4];

		xnos=[xno1,xno2,xno3];
    ynos=[yno1,yno2,yno3];
    znos=[zno1,zno2,zno3];
#     xlabel ('{\it x}')
#     ylabel ('{\it y}')
#     zlabel ('{\it z}')
#     fill3(xnos,ynos,znos,'w');
#     axis equal
#     axis([xmin xmax ymin ymax zmin zmax])
#     hold on;
    xno=sum(xnos)/3;
    yno=sum(ynos)/3;
    zno=sum(znos)/3;
#     text(xno,yno,zno,num2str(i))
    NOS[i,:]=[i,xno,yno,zno];
end;
# view(3)
# hold off;
return NOS
end

function cal_HeG(NOS,NOS_GEO,ELEM,FR,CW,qsi,w,inc)
nelem=size(ELEM,1); # N�mero de elementos (n�mero de linhas da
#  matriz ELEM)

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
		    no4=ELEM[j,5]; # Quarto n� geom�trico do elemento

		    x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
		    y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1
		    z1=NOS_GEO[no1,4]; # Coordenada z do n� geom�trico 1

		    x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
		    y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2
		    z2=NOS_GEO[no2,4]; # Coordenada x do n� geom�trico 2

		    x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 3
		    y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 3
		    z3=NOS_GEO[no3,4]; # Coordenada z do n� geom�trico 3

		    x4=NOS_GEO[no4,2]; # Coordenada x do n� geom�trico 4
		    y4=NOS_GEO[no4,3]; # Coordenada y do n� geom�trico 4
		    z4=NOS_GEO[no4,4]; # Coordenada z do n� geom�trico 4

		    n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio
		    # normal ao elemento
		        if i==j # O ponto fonte pertence ao elemento
		            G[i,j],H[i,j]=calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xd,yd,zd,qsi,w,FR,CW); # Integra��o singular
		        else # O ponto fonte n�o pertence ao elemento
		            G[i,j],H[i,j]=calcula_HeGns(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xd,yd,zd,n,qsi,w,FR,CW); # Integra��o
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

function calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3)
# Function que calcula o vetor unit�rio normal ao elemento

v1 = [x3,y3,z3] - [x2,y2,z2]; # vetor formado pela aresta 32 do elemento
v2 = [x1,y1,z1] - [x2,y2,z2]; # vetor formado pela aresta 12 do elemento
n = cross(v1, v2); # Produto vetorial entre v1 e v2 (vetor normal ao
                           # elemento)
n = n./norm(n); # vetor unit�rio normal ao elemento

return n
end

function calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xd,yd,zd,qsi,w,FR,CW)

# Integra��o singular das matrizes H e G. No caso da matriz G, o elemento
# quadrilateral � decomposto em quatro quadril�teros degenerados na forma
# de tri�ngulos. Os dois primeiros n�s destes quadril�teros s�o
# coincidentes (formam um s� v�rtice do tri�ngulo) e correspondem ao ponto
# onde existe a singularidade, ou seja, ao ponto fonte que se localiza no
# centr�ide do elemento. Isto faz com que haja uma concentra��o de pontos
# de integra��o junto � singularidade, al�m do jacobiano ser igual a zero
# na singularidade. No caso da matriz H, a integra��o � anal�tica e sempre
# ser� igual a -1/2.

npg=length(qsi); # N�mero de pontos de integra��o
g = 0.0; # inicializa��o da matriz G

for kk=1:4
    x3t=xd; # coordenada x do primeiro n� do quadrilatero desgenerado
    y3t=yd; # coordenada y do primeiro n� do quadrilatero desgenerado
    z3t=zd; # coordenada z do primeiro n� do quadrilatero desgenerado
    x4t=xd; # coordenada x do segundo n� do quadrilatero desgenerado
    y4t=yd; # coordenada y do segundo n� do quadrilatero desgenerado
    z4t=zd; # coordenada z do segundo n� do quadrilatero desgenerado

    if(kk==1) # Terceiro e quarto n�s do primeiro quadril�tero desgenerado
        x1t=x1; # coordenada x do terceiro n� do quadrilatero desgenerado
        y1t=y1; # coordenada y do terceiro n� do quadrilatero desgenerado
        z1t=z1; # coordenada z do terceiro n� do quadrilatero desgenerado
        x2t=x2; # coordenada x do quarto n� do quadrilatero desgenerado
        y2t=y2; # coordenada y do quarto n� do quadrilatero desgenerado
        z2t=z2; # coordenada z do quarto n� do quadrilatero desgenerado

    elseif(kk==2) # Terceiro e quarto n�s do segundo quadril�tero
                                  # desgenerado
        x1t=x2; # coordenada x do terceiro n� do quadrilatero desgenerado
        y1t=y2; # coordenada y do terceiro n� do quadrilatero desgenerado
        z1t=z2; # coordenada z do terceiro n� do quadrilatero desgenerado
        x2t=x3; # coordenada x do quarto n� do quadrilatero desgenerado
        y2t=y3; # coordenada y do quarto n� do quadrilatero desgenerado
        z2t=z3; # coordenada z do quarto n� do quadrilatero desgenerado

    elseif(kk==3) # Terceiro e quarto n�s do terceiro quadril�tero
                                  # desgenerado
        x1t=x3; # coordenada x do terceiro n� do quadrilatero desgenerado
        y1t=y3; # coordenada y do terceiro n� do quadrilatero desgenerado
        z1t=z3; # coordenada z do terceiro n� do quadrilatero desgenerado
        x2t=x4; # coordenada x do quarto n� do quadrilatero desgenerado
        y2t=y4; # coordenada y do quarto n� do quadrilatero desgenerado
        z2t=z4; # coordenada z do quarto n� do quadrilatero desgenerado

    elseif(kk==4) # Terceiro e quarto n�s do quarto quadril�tero
                                 # desgenerado
        x1t=x4; # coordenada x do terceiro n� do quadrilatero desgenerado
        y1t=y4; # coordenada y do terceiro n� do quadrilatero desgenerado
        z1t=z4; # coordenada z do terceiro n� do quadrilatero desgenerado
        x2t=x1; # coordenada x do quarto n� do quadrilatero desgenerado
        y2t=y1; # coordenada y do quarto n� do quadrilatero desgenerado
        z2t=z1; # coordenada z do quarto n� do quadrilatero desgenerado
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
            J = calc_jacobiano(x1t,y1t,z1t,x2t,y2t,z2t,x3t,y3t,z3t,x4t,y4t,z4t,qsi[ii],qsi[jj]); # jacobiano(varia ao longo
                         #  do elemento desgenerado)
            Tast,qast = calc_solfund(xd,yd,zd, xc, yc, zc, [0 0 0],FR,CW); # Sol.
            # fudamental de temperatura
            g = g + complex(w[ii] * w[jj] * J,0) * Tast; # integra��o num�rica da matriz G
        end
    end
end

h=complex(1/2,0); # Integra��o anal�tica da matriz H

return g,h
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
function calc_dfforma(qsi, eta)

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

function calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,qsi,eta)
dNdqsi, dNdeta = calc_dfforma(qsi,eta); # Calcula a derivada das fun��es
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

function calc_solfund(x,y,z,xd,yd,zd,n,FR,CW)

rx=x-xd;
ry=y-yd;
rz=z-zd;

r =sqrt(rx^2+ry^2+rz^2);
ZW=complex(0.,-FR*r/CW);
U=exp(ZW)/complex(4*pi*r,0.);
drdn=(rx*n[1] + ry*n[2] + rz*n[3])/r;
Q=(ZW-complex(1.,0.))*U*complex(drdn/r,0.);
return U,Q
end
function calcula_HeGns(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xd,yd,zd,n,qsi,w,FR,CW)
#integra��o n�o singular
n_pint=length(qsi); # N�mero de pontos de integra��o.
g=complex(0,0); # Inicializa o somatorio de g
h=complex(0,0); # Inicializa o somatorio de h


for l=1:n_pint # La�o sobre os pontos de integra��o
    for m=1:n_pint # La�o sobre os pontos de integra��o
        N1,N2,N3,N4=calc_fforma(qsi[l],qsi[m]); #  fun��es de forma
        x=N1*x1+N2*x2+N3*x3+N4*x4; # coordenada x do ponto de integra��o
        y=N1*y1+N2*y2+N3*y3+N4*y4; # coordenada y do ponto de integra��o
        z=N1*z1+N2*z2+N3*z3+N4*z4; # coordenada z do ponto de integra��o
        J = calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,qsi[l],qsi[m]);# jacobiano
        Tast,qast=calc_solfund(x,y,z,xd,yd,zd,n,FR,CW); # Solu��es
        #  fundamentais
        g=g+Tast*complex(w[l]*w[m]*J,0); # Integral da matriz G
        h=h+qast*complex(w[l]*w[m]*J,0); # Integral da matriz H
    end
end
return g,h
end
function aplica_cdc(G,H,CDC)
# Aplica as condi��es de contorno trocando as colunas das matrizes H e G
ncdc = length(CDC[:,1]); # n�mero de linhas da matriz CDC
A=H;
B=G;
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
function monta_Teq(CDC,x)
# Separa fluxo e temperatura

# ncdc = n�mero de linhas da matriz CDC
# T = vetor que cont�m as temperaturas nos n�s
# q = vetor que cont�m o fluxo nos n�s

ncdc = length(CDC[:,1]);
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

function calc_T_pint(PONTOS_int,NOS_GEO,ELEM,T,q,CW,FR,qsi,w,inc)
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
        no4=ELEM[j,5]; # Quarto n� geom�trico do elemento

        x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
        y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1
        z1=NOS_GEO[no1,4]; # Coordenada z do n� geom�trico 1

        x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
        y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2
        z2=NOS_GEO[no2,4]; # Coordenada x do n� geom�trico 2

        x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 3
        y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 3
        z3=NOS_GEO[no3,4]; # Coordenada z do n� geom�trico 3

        x4=NOS_GEO[no4,2]; # Coordenada x do n� geom�trico 4
        y4=NOS_GEO[no4,3]; # Coordenada y do n� geom�trico 4
        z4=NOS_GEO[no4,4]; # Coordenada z do n� geom�trico 4

        n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio
        G_int[i,j],H_int[i,j]=calcula_HeGns(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x_fonte,y_fonte,z_fonte,n,qsi,w,FR,CW); # Chama a functio para calculo de H e G
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
T_pint = - (H_int*T - G_int*q - phi_inc)
#T_pint=-(H_int*T'-G_int*q'-g'); # Vetor que contem a temperatura nos
#      pontos internos
return T_pint
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

function mostra_resultados2(XYZ,tri,T)
nelem=size(tri,1)
x=zeros(3)
y=zeros(3)
z=zeros(3)
zc=zeros(nelem,1)
pc=[zeros(3,3)];
triang=zeros(3,3)
for elem =1:nelem
    no1=tri[elem,2]
    no2=tri[elem,3]
    no3=tri[elem,4]
    x[1]=XYZ[no1,2]
    y[1]=XYZ[no1,3]
    z[1]=XYZ[no1,4]
    x[2]=XYZ[no2,2]
    y[2]=XYZ[no2,3]
    z[2]=XYZ[no2,4]
    x[3]=XYZ[no3,2]
    y[3]=XYZ[no3,3]
    z[3]=XYZ[no3,4]
    triang=[[x[1] y[1] z[1]
    x[2] y[2] z[2]
    x[3] y[3] z[3]]]
    append!(pc,triang)
end
fig = plt.figure()
ax = mp.Axes3D(fig)
q = ar.Poly3DCollection(pc[2:end], linewidths=1)
ax[:add_collection3d](q)
m = cm.ScalarMappable(cmap=cm.jet)
b=m[:to_rgba](T[1:nelem])
q[:set_facecolor](b[:,1:3])
m[:set_array]([minimum(T),maximum(T)])
m[:set_clim](vmin=minimum(T),vmax=maximum(T))
plt.colorbar(m, orientation="vertical",shrink=0.9)
ax[:set_xlabel]("x")
ax[:set_ylabel]("y")
ax[:set_zlabel]("z")
ax[:set_xlim](minimum(XYZ[:,2]),maximum(XYZ[:,2]))
ax[:set_ylim](minimum(XYZ[:,3]),maximum(XYZ[:,3]))
ax[:set_zlim](minimum(XYZ[:,4]),maximum(XYZ[:,4]))
ax[:set_aspect]("equal")
ax[:view_init](elev=18., azim=43.)
return
end
