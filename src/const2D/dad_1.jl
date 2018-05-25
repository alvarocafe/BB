function dad_0(L=1,ne=40)
PONTOS = [1 0 0
	  2 L 0
	  3 L L
	  4 0 L];
SEGMENTOS = [1 1 2 0
	     2 2 3 0
	     3 3 4 0
	     4 4 1 0];
MALHA = [1 ne
	 2 ne
	 3 ne
	 4 ne];
CCSeg = [1 1 0
	 2 0 1
	 3 1 0
	 4 0 0];
return PONTOS, SEGMENTOS, MALHA, CCSeg
end

function dad_1(ne,n_pint,FR)
#Entrada de dados para a comparação com o BEM Isogeometrico
raio = 0.5
PONTOS  = [1   -raio  0 # cilindro com centro em (0.5,0.5) de raio 0.5
          2    raio  0];
#  Segmentos que definem a geometria
#   SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final
#                                                   Raio, tipo do elemento]
#  Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
#                           inicial para o ponto final)
#                    < 0 -> O centro � � direita do segmento (do ponto
#                           inicial para o ponto final)
#                    = 0 -> O segmento � uma linha reta
#  Tipo do elemento = 1 -> Elemento quadrático contínuo
#                   = 2 -> Elemento quadrático descontínuo
#                   = 3 -> Elemento linear contínuo
SEGMENTOS = [1 1 2 -raio  #raio 0.5
             2 2 1 -raio];
# Matriz para defini��o da malha
# MALHA =[numero do segmento, numero de elementos no segmento]
MALHA = [1 ne
         2 ne];
CCSeg=[1 1 1 0
       2 1 1 0];
#FR = 1  #Frequencia analisada
CW = 343  #Velocidade de propagação da onda
#Construindo os pontos internos
#n_pint = 50
PONTOS_int = zeros(n_pint-1,3)
dx = 0
dy = 0.5
passo = 10/(n_pint)
phi_analytical = complex(zeros(size(PONTOS_int,1),1));
for i = 1:n_pint-1  # Para n_pint pontos internos
  PONTOS_int[i,:] = [i  dx dy+i*passo]
  bh1 = SpecialFunctions.besselh(1,2,(FR/CW)*raio);
  bh0 = SpecialFunctions.besselh(0,2,(FR/CW)*PONTOS_int[i,3]);
  phi_analytical[i,1] = (CW/FR).*(bh0./bh1);		#solucao anali­tica pela separacao de variaveis em coordenadas cilindricas da equacao de Helmholtz
end
#Incluimos as fontes concentradas e as ondas incidentes
fc = [0 3*0.5/2 -0.5/2 1];
finc = [0 -1 0 1];
return PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc,phi_analytical
end

function dad_vibcyl(ne,k,n_pint=50,raio=0.5,delta=10)
#Entrada de dados para a comparação com o BEM Isogeometrico
#raio = 0.5
PONTOS  = [1   -raio  0 # cilindro com centro em (0.5,0.5) de raio 0.5
          2    raio  0];
#  Segmentos que definem a geometria
#   SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final
#                                                   Raio, tipo do elemento]
#  Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
#                           inicial para o ponto final)
#                    < 0 -> O centro � � direita do segmento (do ponto
#                           inicial para o ponto final)
#                    = 0 -> O segmento � uma linha reta
#  Tipo do elemento = 1 -> Elemento quadrático contínuo
#                   = 2 -> Elemento quadrático descontínuo
#                   = 3 -> Elemento linear contínuo
SEGMENTOS = [1 1 2 -raio  #raio 0.5
             2 2 1 -raio];
# Matriz para defini��o da malha
# MALHA =[numero do segmento, numero de elementos no segmento]
MALHA = [1 ne
         2 ne];
CCSeg=[1 1 1 0
       2 1 1 0];

# Create the domain points
passo = delta/(n_pint - 1);
iter = 1;
for i = 2:n_pint
	for j = 2:n_pint
		if ((i*passo)^2 + (j*passo)^2) > 30*raio
			iter +=4
		end
	end
end
PONTOS_int = zeros(iter-1,3)
phi_analytical = complex(zeros(iter-1,1))
iter = 1;
for i = 2:n_pint
	for j = 2:n_pint
		if ((i*passo)^2 + (j*passo)^2) > 30*raio
			PONTOS_int[iter,:] = [iter  (i-1)*passo (j-1)*passo]
			  bh1 = SpecialFunctions.besselh(1,2,(FR/CW)*raio);
			  bh0 = SpecialFunctions.besselh(0,2,(FR/CW)*sqrt((PONTOS_int[iter,3])^2 + (PONTOS_int[iter,2]^2)));
			  phi_analytical[iter:iter+3,1] = (CW/FR).*(bh0./bh1);		#solucao anali­tica pela separacao de variaveis em coordenadas cilindricas da equacao de Helmholtz
			PONTOS_int[iter+1,:] = [iter+1  (i-1)*passo -(j-1)*passo]
			PONTOS_int[iter+2,:] = [iter+2  -(i-1)*passo -(j-1)*passo]
			PONTOS_int[iter+3,:] = [iter+3  -(i-1)*passo (j-1)*passo]
			iter +=4
		end
	end
end
#Incluimos as fontes concentradas e as ondas incidentes
fc = [0 3*0.5/2 -0.5/2 1];
finc = [0 -1 0 1];
return PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int,fc,finc,phi_analytical
end


function dad_lev_res(L=10,ne = 2,n_pint=100,delta=100)
# This test case consists of a rectangular reflector and an acoustic source so that acoustic levitation is possible on the modal nodes between the source and the resonator.

PONTOS  = [1 0 0 ;
  	   2 L 0 ;
	   3 L L/10 ;
	   4 0 L/10 ];
#  Segmentos que definem a geometria
#   SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final
#                                                   Raio, tipo do elemento]
#  Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
#                           inicial para o ponto final)
#                    < 0 -> O centro � � direita do segmento (do ponto
#                           inicial para o ponto final)
#                    = 0 -> O segmento � uma linha reta
#  Tipo do elemento = 1 -> Elemento quadrático contínuo
#                   = 2 -> Elemento quadrático descontínuo
#                   = 3 -> Elemento linear contínuo
SEGMENTOS = [1 1 4 0
             2 4 3 0
	     3 3 2 0
	     4 2 1 0];
# Matriz para defini��o da malha
# MALHA =[numero do segmento, numero de elementos no segmento]
MALHA = [1 ne
         2 ne
	 3 ne
	 4 ne];
CCSeg=[ 1 1 1 1
	2 1 1 0
	3 1 1 0
	4 1 1 0]	

CW = 344*100  #Velocidade de propagação da onda [cm/s]
FR = 40000  #Frequencia analisada corresponde a um comprimento de onda de 1 cm
cp = CW/FR	# Wavelength
k = 1/cp		# Wavenumber

# Adding concentrated sources and incident waves
# For these terms, the data structure is:
	# fc = [Q x y z]; where Q is the strength of the concentrated source and x,y,z are the position of the source  
	# finc = [A x y z]; where A is the amplitude of the plane wave and x,y,z are the direction d = (x,y,z) of propagation such that |d| = 1 
fc = [0 3*0.5/2 -0.5/2 1];
d = 5;
finc = [0 5 -d 0];
phi_analytical = []

# External points creation
passo = delta/(n_pint - 1);
iter = 1;
for i = 2:n_pint
	for j = 2:n_pint
		if ((i*passo)^2 + (j*passo)^2) > 1.5*100^2
			iter +=4
		end
	end
end
PONTOS_int = zeros(iter-1,3)
iter = 1;
for i = 2:n_pint
	for j = 2:n_pint
		if ((i*passo)^2 + (j*passo)^2) > 1.5*100^2
			PONTOS_int[iter,:] = [iter  (i-1)*passo (j-1)*passo]
			PONTOS_int[iter+1,:] = [iter+1  (i-1)*passo -(j-1)*passo]
			PONTOS_int[iter+2,:] = [iter+2  -(i-1)*passo -(j-1)*passo]
			PONTOS_int[iter+3,:] = [iter+3  -(i-1)*passo (j-1)*passo]
			iter +=4
		end
	end
end
return PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc,phi_analytical,k
end

function dad_rocket(ne = 2,delta=1000,n_pint=100)
# This test case consists of a rectangular reflector and an acoustic source so that acoustic levitation is possible on the modal nodes between the source and the resonator.

PONTOS  = [1 0 0 ;
  	   2 10 0 ;
	   3 10 100 ;
	   4 0 100 ];
#  Segmentos que definem a geometria
#   SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final
#                                                   Raio, tipo do elemento]
#  Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
#                           inicial para o ponto final)
#                    < 0 -> O centro � � direita do segmento (do ponto
#                           inicial para o ponto final)
#                    = 0 -> O segmento � uma linha reta
#  Tipo do elemento = 1 -> Elemento quadrático contínuo
#                   = 2 -> Elemento quadrático descontínuo
#                   = 3 -> Elemento linear contínuo
SEGMENTOS = [1 1 4 0
             2 4 3 0
	     3 3 2 0
	     4 2 1 0];
# Matriz para defini��o da malha
# MALHA =[numero do segmento, numero de elementos no segmento]
MALHA = [1 ne
         2 ne
	 3 ne
	 4 ne];
CCSeg=[ 1 1 1 0
	2 1 1 0
	3 1 1 0
	4 1 0 1]	

CW = 344  #Velocidade de propagação da onda [m/s]
FR = 2*pi*20

# Adding concentrated sources and incident waves
# For these terms, the data structure is:
	# fc = [Q x y z]; where Q is the strength of the concentrated source and x,y,z are the position of the source  
	# finc = [A x y z]; where A is the amplitude of the plane wave and x,y,z are the direction d = (x,y,z) of propagation such that |d| = 1 
fc = [0 3*0.5/2 -0.5/2 1];
d = 5;
finc = [1 5 -d 0];
#phi_analytical = []

# External points creation
#Construindo os pontos internos
#delta = 1000 # Max distance from domain point
#n_pint =100 # Number of internal points
passo = delta/(n_pint - 1);
iter = 1;
for i = 2:n_pint
	for j = 2:n_pint
		if ((i*passo)^2 + (j*passo)^2) > 1.5*100^2
			iter +=4
		end
	end
end
PONTOS_int = zeros(iter-1,3)
iter = 1;
for i = 2:n_pint
	for j = 2:n_pint
		if ((i*passo)^2 + (j*passo)^2) > 1.5*100^2
			PONTOS_int[iter,:] = [iter  (i-1)*passo (j-1)*passo]
			PONTOS_int[iter+1,:] = [iter+1  (i-1)*passo -(j-1)*passo]
			PONTOS_int[iter+2,:] = [iter+2  -(i-1)*passo -(j-1)*passo]
			PONTOS_int[iter+3,:] = [iter+3  -(i-1)*passo (j-1)*passo]
			iter +=4
		end
	end
end
# Adding concentrated sources and incident waves
# For these terms, the data structure is:
	# fc = [Q x y z]; where Q is the strength of the concentrated source and x,y,z are the position of the source  
	# finc = [A x y z]; where A is the amplitude of the plane wave and x,y,z are the direction d = (x,y,z) of propagation such that |d| = 1 
fc = [0 3*0.5/2 -0.5/2 1];
d = 5;
finc = [0 5 -d 0];

return PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc
end
