function dad_const(p1,p2)
  #Essa é a função para produzir um elemento linear constante
  #  Segmentos que definem a geometria
  #   SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final
  #                                                   Raio, tipo do elemento]
  #  Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
  #                           inicial para o ponto final)
  #                    < 0 -> O centro � � direita do segmento (do ponto
  #                           inicial para o ponto final)
  #                    = 0 -> O segmento � uma linha reta
  PONTOS = [1 p1[1] p1[2]
            2 p2[1] p2[2]]
  SEGMENTOS = [1 1 2 0];  #Raio igual a zero, elemento retilíneo

  # Matriz para definição da malha

  # MALHA =[numero do segmento, numero de elementos no segmento]

  ne=1;
  MALHA = [1 ne];
  CCSeg=[1 1 0];

  #Geramos as matrizes de pontos geometricos, nós, elementos e condições de contorno
  NOS_GEO = [1 p1[1] p1[2]
             2 p2[1] p2[2]]
  NOS = [1 (p2[1]-p1[1])/2 (p2[2]-p1[2])/2]
  ELEM = [1 1 2]
  CDC = [1 0 1]

  return NOS_GEO,NOS,ELEM,CDC
end

function dad_quad(p1,p2,p3)
  #Essa função cria apenas um elemento quadrático 2D para testes
  NOS = [1  p1[1] p1[2]
         2  p2[1] p2[2]
         3  p3[1] p3[2]]
  ELEM = [1 1 2 3]
  CDC = [1 1 1]
  return NOS, ELEM, CDC
end

function dad_iso(ne)
  #Essa é a função para produzir um elemento NURBS
  #  Segmentos que definem a geometria
  #   SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final
  #                                                   Raio, tipo do elemento]
  #  Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
  #                           inicial para o ponto final)
  #                    < 0 -> O centro � � direita do segmento (do ponto
  #                           inicial para o ponto final)
  #                    = 0 -> O segmento � uma linha reta
  PONTOS = [1 p1[1] p1[2]
            2 p2[1] p2[2]]
  # SEGMENTOS = [1 1 2 0];  #Raio igual a zero, elemento retilíneo
  SEGMENTOS = [1 1 2 0.5];  #Raio = 0.5, elemento curvo circular
  # Matriz para definição da malha
  # MALHA =[numero do segmento, numero de elementos no segmento]
  ne=1;
  MALHA = [1 ne];
  CCSeg=[1 1 1];
  # Gerando a curva NURBS
  crv = format_dad_iso(PONTOS,SEGMENTOS,MALHA)
  dcrv=map(x->nrbderiv(x),crv)
  n = length(crv);	# N�mero total de elementos
#  p=1;#refinamento p
#  for i=1:n
#      degree=crv[i].order-1
#      #	println(crv[i].knots)
#      #	println(crv[i].coefs)
#      coefs,knots = bspdegelev(degree,crv[i].coefs,crv[i].knots,p)
#      #	println(knots)
#      #	println(coefs)
#      crv[i] = nrbmak(coefs,knots)
#  end
#  h=10;#refinamento h
#  for i=1:n
#    novosnos=linspace(0,1,h+2)
#    degree=crv[i].order-1
#    coefs,knots = bspkntins(degree,crv[i].coefs,crv[i].knots,novosnos[2:end-1])
#    crv[i] = nrbmak(coefs,knots)
#  end
  z=0;
  for k=1:n
      for i=1:crv[k].number
          z=z+1
      end
  end
  numcurva=zeros(Integer,z)
  collocPts=zeros(z)
  CDC=zeros(z,3)
  collocCoord=zeros(z,3)
  z=0;
  nnos=zeros(Integer,n)
  for k=1:n
      p=crv[k].order-1;
      nnos[k]=crv[k].number;
    valorCDC=CCSeg[k,3];
    tipoCDC=CCSeg[k,2];
    for i=1:crv[k].number
        z=z+1;
        numcurva[z]=k;
        collocPts[z]=sum(crv[k].knots[(i+1):(i+p)])/p;
        if(i==2)
            collocPts[z-1]=(collocPts[z]+collocPts[z-1])/2;
        end
        if(i==nnos[k])
            collocPts[z]=(collocPts[z]+collocPts[z-1])/2;
        end

       CDC[z,:] = [z,tipoCDC,valorCDC];
    end
  end
  nnos2=cumsum([0 nnos'],2);

  E=zeros(length(collocPts),length(collocPts));
  for i=1:length(collocPts)
    collocCoord[i,:]=nrbeval(crv[numcurva[i]], collocPts[i]);
    B, id = nrbbasisfun(crv[numcurva[i]],collocPts[i])
    E[i,id+nnos2[numcurva[i]]]=B
  end
    # plot(collocCoord[:,1],collocCoord[:,2])
    # legend('Curva resultante','Polígono de controle','Pontos de controle','Pontos fonte')
  return collocCoord,nnos2,crv,dcrv,CDC,E
end

function dad_0(ne)
    #Entrada de dados para testar o programa de elementos constantes
    PONTOS  = [1   0  0
               2   1  0
               3   1  1
               4   0  1];
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
    SEGMENTOS = [1 1 2  0  #raio 0
                 2 2 3  0
                 3 3 4  0
                 4 4 1  0];
    # Matriz para defini��o da malha
    # MALHA =[numero do segmento, numero de elementos no segmento]
    MALHA = [1 ne
             2 ne
             3 ne
             4 ne];
    CCSeg=[1 1 0
           2 0 1
           3 1 0
           4 1 0];
    AFR = 1;
    CW = 1;
    fc = [1 3*0.5/2 -0.5/2 1];
    finc = [1 -1 0 1];
    return PONTOS, SEGMENTOS, MALHA, CCSeg, AFR, CW,fc,finc
end

function dad_1(ne,FR)
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
n_pint = 50
PONTOS_int = zeros(n_pint,3)
dx = 1.0
dy = 0
passo = 0.05
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

function dad_quad_cyl(ne)
#Esse programa gera a matriz NOS, ELEM e CDC para o elemento quadratico
#para um cilindro vibrante com ne elementos para discretizar o círculo
# Esse programa gera o cilindro com o centro em (0,0). É possível alterar pra outros centros
raio = 0.5
nnos = 2*ne
no = 1
ELEM = zeros(ne,4)
CDC = zeros(ne,4)
NOS = zeros(2*ne,3)
for i = 1:ne
  if i == ne
    ELEM[i,:] = [i no no+1 1]
    NOS[no:no+1,:] = [no raio*cos((2*pi/(nnos))*(no-1)) -raio*sin((2*pi/(nnos))*(no-1))
                      no+1 raio*cos((2*pi/(nnos))*(no)) -raio*sin((2*pi/(nnos))*(no))]
  else
    ELEM[i,:] = [i no no+1 no+2]
    NOS[no:no+2,:] = [no raio*cos((2*pi/(nnos))*(no-1)) -raio*sin((2*pi/(nnos))*(no-1))
                      no+1 raio*cos((2*pi/(nnos))*(no)) -raio*sin((2*pi/(nnos))*(no))
                      no+2 raio*cos((2*pi/(nnos))*(no+1)) -raio*sin((2*pi/(nnos))*(no+1))]
    no = no+2
  end
end
CDC = ones(ne,4)
CDC[:,4] = 0
FR = 1  #Frequencia analisada
CW = 1  #Velocidade de propagação da onda
#Construindo os pontos internos
n_pint = 50
PONTOS_int = zeros(n_pint,3)
dx = 1.0
dy = 0
passo = 0.05
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

  return NOS, NOS, ELEM, CDC, PONTOS_int, phi_analytical, CW, FR, fc, finc
end

function dad_ext(n_el)
n_el = n_el/2
#Entrada de dados para problema do cilindro pulsante BEM
#Autor: Alvaro Campos Ferreira
c = [0 0]; #Centro do furo [x,y]
raio = 0.5; #Raio em metros
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

# CW=340; # Velocidade de propagaÃ§Ã£o de onda no meio (ar = 344 m/s)
# wb=pi*CW/10; # Primeira frequÃªncia natural do problema (soluÃ§Ã£o analÃ­tica)
CW = 1
wb = 1
freq_ini=0.1*wb; # FrequÃªncia inicial
freq_final=1.1*wb; # FrequÃªncia final
nfreq=10; # NÃºmero de frequÃªncias que serÃ£o analisadas
AFR=linspace(freq_ini,freq_final,nfreq); # FrequÃªncias que serÃ£o
                             #    analisadas (em rad/s)
AFR=wb;

#Rotina para carga concentrada
fc=[0 0 0 0];
finc = [0 0 0 0];
#CriaÃ§Ã£o dos pontos externos
flag = 1; #Esse flag serÃ¡ usado na visualizaÃ§Ã£o de dados, indica que os pontos estÃ£o em uma linha reta em x ou y
PONTOS_int = zeros(10,3);
xint = 0.51;
passo = 0.1;
phi_analytical = complex(zeros(size(PONTOS_int,1),length(AFR)));
for i = 1:size(PONTOS_int,1)
	PONTOS_int[i,1] = i;
	PONTOS_int[i,2] = xint;
	PONTOS_int[i,3] = 0;
	#SoluÃ§Ã£o analÃ­tica para o cilindro pulsante infinito
	for j = 1:length(AFR)
		k = AFR[j]./CW;
    bh1 = SpecialFunctions.besselh(1,2,k.*raio);
		bh0 = SpecialFunctions.besselh(0,2,k.*xint);
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
return PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, AFR, CW, fc, finc,phi_analytical
end

function dad_triconst(p1,p2,p3)
  NOS_GEO = [1 p1[1] p1[2]  p1[3]
             2 p2[1] p2[2]  p2[3]
             3 p3[1] p3[2]  p3[3]]
  NOS = [1 (p1[1] + p2[1] + p3[1])/3 (p1[2] + p2[2] + p3[2])/3  (p1[3] + p2[3] + p3[3])/3]
  ELEM = [1 1 2 3]
  CDC = [1 1 1]
  return NOS_GEO,NOS,ELEM,CDC
end

function dad_quadlin(p1,p2,p3,p4)
    NOS_GEO = [1 p1[1] p1[2]  p1[3]
               2 p2[1] p2[2]  p2[3]
               3 p3[1] p3[2]  p3[3]
               4 p4[1] p4[2]  p4[3]]
    NOS = [1 (p1[1] + p2[1] + p3[1]+ p4[1])/4 (p1[2] + p2[2] + p3[2]+ p4[2])/4  (p1[3] + p2[3] + p3[3]+ p4[3])/4]
    ELEM = [1 1 2 3 4]
    CDC = [1 1 1]
    return  NOS_GEO,NOS,ELEM,CDC
end

function dad_triquad(p1,p2,p3)
  NOS_GEO = [1 p1[1] p1[2]  p1[3]
             2 p2[1] p2[2]  p2[3]
             3 p3[1] p3[2]  p3[3]
             4 (p1[1]+p2[1])/2 (p1[2]+p2[2])/2  (p1[3]+p2[3])/2
             5 (p2[1]+p3[1])/2 (p2[2]+p3[2])/2  (p2[3]+p3[3])/2
             6 (p1[1]+p3[1])/2 (p1[2]+p3[2])/2  (p1[3]+p3[3])/2]
  NOS = [1 (p1[1] + p2[1] + p3[1])/3 (p1[2] + p2[2] + p3[2])/3  (p1[3] + p2[3] + p3[3])/3]
  ELEM = [1 1 4 2 5 3 6]
  CDC = [1 1 1]
  return NOS_GEO,NOS,ELEM,CDC
end

function dad_quadquad(p1,p2,p3,p4)
  NOS_GEO = [1 p1[1] p1[2]  p1[3]
             2 p2[1] p2[2]  p2[3]
             3 p3[1] p3[2]  p3[3]
             4 p4[1] p4[2]  p4[3]
             5 (p1[1]+p2[1])/2 (p1[2]+p2[2])/2  (p1[3]+p2[3])/2
             6 (p2[1]+p3[1])/2 (p2[2]+p3[2])/2  (p2[3]+p3[3])/2
             7 (p3[1]+p4[1])/2 (p3[2]+p4[2])/2  (p3[3]+p4[3])/2
             8 (p4[1]+p1[1])/2 (p4[2]+p1[2])/2  (p4[3]+p1[3])/2]

  NOS = [1 (p1[1] + p2[1] + p3[1] + p4[1])/4 (p1[2] + p2[2] + p3[2] + p4[2])/4  (p1[3] + p2[3] + p3[3] + p4[3])/4]
  ELEM = [1 1 5 2 6 3 7 4 8]
  CDC = [1 1 1]
    return NOS_GEO,NOS,ELEM,CDC
end
