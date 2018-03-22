function dad_1()
#Entrada de dados para a comparação com o BEM Isogeometrico
PONTOS  = [1   0  0 # cilindro com centro em (0.5,0.5) de raio 0.5
          2    1  0];

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

et = 1;
SEGMENTOS = [1 1 2  .5  #raio 0.5
         2 2 1  .5];

# Matriz para defini��o da malha

# MALHA =[numero do segmento, numero de elementos no segmento]

ne=25;
MALHA = [1 ne
         2 ne];
CCSeg=[1 1 0
    2 0 1];

 k=1;
NPX = 1;
NPY = 1;
AFR = 1;
CW = 1;
PONTOS_int = [1 0.5 0.5];
fc = [1 3*0.5/2 -0.5/2 1];
finc = [1 -1 0 1];
return PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, AFR, CW,fc,finc
end


function dad_ext(n_el=100)
#Entrada de dados para problema do cilindro pulsante BEM
#Autor: Alvaro Campos Ferreira
c = [0 0]; #Centro do furo [x,y]
raio = 0.1; #Raio em metros
ne_furos = n_el; #Numero de elementos por quarto de furo

PONTOS=	[1	c[1,1]-raio	c[1,2] 	#Coordenadas dos pontos fÃ­sicos
       	2	c[1,1]		c[1,2]+raio
       	3	c[1,1]+raio	c[1,2]
       	4	c[1,1]		c[1,2]-raio];

SEGMENTOS=	[1 1 2 -raio		#Segmentos de reta com curvatura de -raio
          	 2 2 3 -raio
                 3 3 4 -raio
            	 4 4 1 -raio];
# CondiÃ§Ãµes de contorno nos segmentos
# CCSeg=[no do segmento, tipo da CDC, valor REAL da CDC, valor IMAG da CDC]
# tipo da CDC = 0 => a pressÃ£o Ã© conhecida
# tipo da CDC = 1 => O fluxo Ã© conhecido
CCSeg=[1 1 1 0			#CondiÃ§Ãµes de contorno para cada segmento
       2 1 1 0
       3 1 1 0
       4 1 1 0];
MALHA=[1 ne_furos		#NÃºmero de elementos por segmento
       2 ne_furos
       3 ne_furos
       4 ne_furos];

CW=340; # Velocidade de propagaÃ§Ã£o de onda no meio (ar = 344 m/s)
wb=pi*CW/10; # Primeira frequÃªncia natural do problema (soluÃ§Ã£o analÃ­tica)
freq_ini=0.1*wb; # FrequÃªncia inicial
freq_final=1.1*wb; # FrequÃªncia final
nfreq=10; # NÃºmero de frequÃªncias que serÃ£o analisadas
AFR=linspace(freq_ini,freq_final,nfreq); # FrequÃªncias que serÃ£o
                             #    analisadas (em rad/s)
#  AFR=wb;

#Rotina para carga concentrada
fc=[0 0 0 0];
finc = [0 0 0 0];
#CriaÃ§Ã£o dos pontos externos
flag = 1; #Esse flag serÃ¡ usado na visualizaÃ§Ã£o de dados, indica que os pontos estÃ£o em uma linha reta em x ou y
PONTOS_int = zeros(10,3);
xint = 0.11;
passo = 0.1;
phi_analytical = complex(zeros(size(PONTOS_int,1),length(AFR)));
for i = 1:size(PONTOS_int,1)
	PONTOS_int[i,1] = i;
	PONTOS_int[i,2] = xint;
	PONTOS_int[i,3] = 0;
	#SoluÃ§Ã£o analÃ­tica para o cilindro pulsante infinito
	for j = 1:length(AFR)
		k = AFR[j]./CW;
		bh0 = SpecialFunctions.besselh(0,2,k.*xint);
		bh1 = SpecialFunctions.besselh(1,2,k.*raio);
		phi_analytical[i,j] = (1./k).*(bh0./bh1);		#soluÃ§Ã£o analÃ­tica pela separaÃ§Ã£o de variÃ¡veis em coordenadas cilindricas da equaÃ§Ã£o de Helmholtz
	end
	xint = xint + passo;	#incremento de x

end

#passo_theta = 2*pi/10;
#theta[:,1] = 0:passo_theta:10*passo_theta;
#break
#R[:,1] = PONTOS_int[:,2];
#for i=2:size(theta)
#	R[:,i] = PONTOS_int[:,2];
#	theta[:,i-1] = 0:passo_theta:10*passo_theta;
#end
#size(R)
#size(theta')
#[X,Y]=pol2cart(theta',R);
return PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, AFR, CW, fc, finc
end
