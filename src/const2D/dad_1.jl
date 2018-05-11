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
PONTOS_int = zeros(n_pint,3)
dx = 0.5
dy = 0
passo = 1/(n_pint)
phi_analytical = complex(zeros(size(PONTOS_int,1),1));
for i = 1:n_pint  # Para n_pint pontos internos
  PONTOS_int[i,:] = [i  dx+i*passo dy]
  bh1 = SpecialFunctions.besselh(1,2,(FR/CW)*raio);
  bh0 = SpecialFunctions.besselh(0,2,(FR/CW)*PONTOS_int[i,2]);
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
#FR = 1  #Frequencia analisada
#CW = 343  #Velocidade de propagação da onda
#Construindo os pontos internos
passo = delta/(n_pint - 1);
iter = 1;
for i = 2:n_pint
	for j = 2:n_pint
		if ((i*passo)^2 + (j*passo)^2) > 1.5*raio^2
			iter +=4
		end
	end
end
println(iter)
PONTOS_int = zeros(iter-1,3)
phi_analytical = complex(zeros(iter-1))
bh1 = SpecialFunctions.besselh(1,2,k*raio);
iter = 1;
for i = 2:n_pint
	for j = 2:n_pint
		if ((i*passo)^2 + (j*passo)^2) > 1.5*raio^2
			PONTOS_int[iter,:] = [iter  (i-1)*passo (j-1)*passo]
			PONTOS_int[iter+1,:] = [iter+1  (i-1)*passo -(j-1)*passo]
			PONTOS_int[iter+2,:] = [iter+2  -(i-1)*passo -(j-1)*passo]
			PONTOS_int[iter+3,:] = [iter+3  -(i-1)*passo (j-1)*passo]
			bh0 = SpecialFunctions.besselh(0,2,k*sqrt((i*passo)^2 + (j*passo)^2));
			phi_analytical[iter] = (1/k).*(bh0./bh1);		#solucao anali­tica pela 
			phi_analytical[iter+1] = phi_analytical[iter];		#solucao anali­tica pela 
			phi_analytical[iter+2] = phi_analytical[iter];		#solucao anali­tica pela 
			phi_analytical[iter+3] = phi_analytical[iter];		#solucao anali­tica pela 
			iter +=4
		end
	end
end

#Incluimos as fontes concentradas e as ondas incidentes
fc = [0 3*0.5/2 -0.5/2 1];
finc = [0 -1 0 1];
return PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int,fc,finc,phi_analytical
end


function dad_lev_res(ne = 2)
# This test case consists of a rectangular reflector and an acoustic source so that acoustic levitation is possible on the modal nodes between the source and the resonator.

PONTOS  = [1 0 0 ;
  	   2 10 0 ;
	   3 10 1 ;
	   4 0 1 ];
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
finc = [1 5 -d 0];
phi_analytical = []

# External points creation


PONTOS_int = []

return PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc,phi_analytical
end

