function dad_0(L=1)
PONTOS = [1 0 0
	  2 L 0
	  3 L L
	  4 0 L];
SEGMENTOS = [1 1 2 0
	     2 2 3 0
	     3 3 4 0
	     4 4 1 0];
ne = 4;
MALHA = [1 ne
	 2 ne
	 3 ne
	 4 ne];
CCSeg = [1 1 0
	 2 0 1
	 3 1 0
	 4 0 0];
crv = format_dad_iso(PONTOS,SEGMENTOS,MALHA)
dcrv=map(x->nrbderiv(x),crv)
n = length(crv);	# N�mero total de elementos
  p=1;#refinamento p
  for i=1:n
      degree=crv[i].order-1
      #	println(crv[i].knots)
      #	println(crv[i].coefs)
      coefs,knots = bspdegelev(degree,crv[i].coefs,crv[i].knots,p)
      #	println(knots)
      #	println(coefs)
      crv[i] = nrbmak(coefs,knots)
  end
  h=10;#refinamento h
  for i=1:n
    novosnos=linspace(0,1,h+2)
    degree=crv[i].order-1
    coefs,knots = bspkntins(degree,crv[i].coefs,crv[i].knots,novosnos[2:end-1])
    crv[i] = nrbmak(coefs,knots)
  end
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


function dad_iso(raio)
  #Essa é a função para produzir um elemento NURBS
  #  Segmentos que definem a geometria
  #   SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final
  #                                                   Raio, tipo do elemento]
  #  Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
  #                           inicial para o ponto final)
  #                    < 0 -> O centro � � direita do segmento (do ponto
  #                           inicial para o ponto final)
  #                    = 0 -> O segmento � uma linha reta
  PONTOS = [1 0 0
            2 1 0]
  # SEGMENTOS = [1 1 2 0];  #Raio igual a zero, elemento retilíneo
  SEGMENTOS = [1 1 2 raio
	       2 2 1 raio];  #Raio = 0.5, elemento curvo circular
  # Matriz para definição da malha
  # MALHA =[numero do segmento, numero de elementos no segmento]
  ne=1;
  MALHA = [1 ne
	   2 ne];
  CCSeg=[1 1 1
	 2 1 1];
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

function dad_helm(L=6)
PONTOS = [1 0 0
	  2 L 0
	  3 L L
	  4 0 L];
SEGMENTOS = [1 1 2 0
	     2 2 3 0
	     3 3 4 0
	     4 4 1 0];
ne = 4;
MALHA = [1 ne
	 2 ne
	 3 ne
	 4 ne];
CCSeg = [1 0 0
	 2 1 0
	 3 1 100
	 4 0 0];
  # Gerando a curva NURBS
  crv = format_dad_iso(PONTOS,SEGMENTOS,MALHA)
  dcrv=map(x->nrbderiv(x),crv)
  n = length(crv);	# N�mero total de elementos
  p=1;#refinamento p
  for i=1:n
      degree=crv[i].order-1
      #	println(crv[i].knots)
      #	println(crv[i].coefs)
      coefs,knots = bspdegelev(degree,crv[i].coefs,crv[i].knots,p)
      #	println(knots)
      #	println(coefs)
      crv[i] = nrbmak(coefs,knots)
  end
  h=10;#refinamento h
  for i=1:n
    novosnos=linspace(0,1,h+2)
    degree=crv[i].order-1
    coefs,knots = bspkntins(degree,crv[i].coefs,crv[i].knots,novosnos[2:end-1])
    crv[i] = nrbmak(coefs,knots)
  end
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
