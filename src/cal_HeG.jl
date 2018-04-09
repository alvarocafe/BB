function cal_HeG(NOS,NOS_GEO,ELEM,CW,FR,fc,qsi,w)
  nelem::Int64=size(ELEM)[1]; # N�mero de elementos (n�mero de linhas da matriz ELEM)
  nnos::Int64=nelem; # N�mero de n�s
  G=complex(zeros(nnos,nnos));
  H=complex(zeros(nnos,nnos));
  q=zeros(nnos,1);  #Influencia de fontes concentradas
  inc = zeros(nnos,1);  #Influencia de ondas incidentes
  qsitelles,Jtelles = telles(qsi,0); # Calculando os pontos de Telles e o Jacobiano

  for i=1:nnos # Laco sobre os pontos fontes
    xd=NOS[i,2]; # Coordenada x do ponto fonte
    yd=NOS[i,3]; # Coordenada y do ponto fonte
    for j=1:nelem # Laco sobre os elementos
      noi::Int64=ELEM[j,2]; # Ponto inicial do elemento
      nof::Int64=ELEM[j,3]; # Ponto final do elemento
      x1=NOS_GEO[noi,2]; # Coordenada x do ponto inicial do elemento
      x2=NOS_GEO[nof,2]; # Coordenada x do ponto final do elemento
      y1=NOS_GEO[noi,3]; # Coordenada y do ponto inicial do elemento
      y2=NOS_GEO[nof,3];  # Coordenada y do ponto final do elemento
      if i==j # O ponto fonte pertence ao elemento
        # g,h = calcula_HeGs(x1,y1,x2,y2,CW,FR);
        g, h = calcula_HeGns(x1,y1,x2,y2,xd,yd,CW,qsitelles,w.*Jtelles,FR);
        #h = h + 0.5
        # println("Diferença entre g e gtelles = ", abs(g-gtelles))
        # println("Diferença entre h e htelles = ", abs(h-htelles))
      else # O ponto fonte n�o pertence ao elemento
        g,h = calcula_HeGns(x1,y1,x2,y2,xd,yd,CW,qsi,w,FR);
      end
      G[i,j] = g
      H[i,j] = h
    end
  end
  # Calculo dos termos da diagonal da matriz H [consideração de corpo a
# temperatura constante]
for m = 1 : nnos
    H[m,m] = 0; #zera a diagonal principal
    for n = 1 : nnos
        if n != m
            H[m,m] = H[m,m] - H[m,n];
        end;
    end;
end;

  return G,H
end

function calcula_HeGns(x1,y1,x2,y2,xd,yd,CW,qsi,w,FR)
  #integração não singular
  n_pint=length(qsi); # Número de pontos de integração.
  g=complex(0,0); # Inicializa o somatorio de g
  h=complex(0,0); # Inicializa o somatorio de h
  L=sqrt((x2-x1)^2+(y2-y1)^2); # Comprimento do elemento
  dgamadqsi=L/2; # Jacobiano

  sx=(x2-x1)/L; # Componente x do vetor tangente
  sy=(y2-y1)/L; # Componente y do vetor tangente
  nx=sy; # Componente x do vetor normal
  ny=-sx; # Componente y do vetor normal

  for kk=1:n_pint # Laço sobre os pontos de integração
    N1,N2 =calc_fforma(qsi[kk]); # Calcula as funções de forma
    x=N1*x1+N2*x2; # Calcula a coordenada x do ponto de integração
    y=N1*y1+N2*y2; # Calcula a coordenada y do ponto de integração

    Tast,qast =calc_solfund(x,y,xd,yd,nx,ny,CW,FR); # Obtemos as soluções fundamentais

    h=h+qast*w[kk]*dgamadqsi; # Integral da matriz H
    g=g+Tast*w[kk]*dgamadqsi; # Integral da matriz G
  end
  return g,h
end

function calcula_HeGs(X1,Y1,X2,Y2,CW,FR)
  # integração singular
  #  THIS SUBROUTINE COMPUTES THE VALUES OF THE DIAGONAL
  #  COEFFICIENTS OF THE G MATRIX

  AX=(X2-X1)/2.;
  AY=(Y2-Y1)/2.;
  SR=sqrt(AX^2+AY^2);
  ZR=real(FR*SR/CW);
  Z=complex(0.,ZR);
  S0C=INBESS(Z);
  G=(SR/Z)*S0C/(pi);
  H=0.5;
  return G, H
end

function INBESS(Z)
  #
  #  THIS SUBROUTINE COMPUTES THE INTEGRAL OF THE MODIFIED BESSEL FUNCTION
  #  OF ZERO ORDER BETWEEN 0 AND Z USING A POLINOMIUM.
  #
  #  THE COEFFICIENTS OF THE POLINOMIUM HAVE BEEN OBTAINED FROM THE
  #  TERMS OF THE SERIES EXPANSION

  S0=0.42278434;
  S1=[0.90734120e+01, 0.72756425E+02, 0.25901019E+03,  0.52398212e+03, 0.68598066e+03, 0.62976149e+03,  0.42826262E+03, 0.22451590E+03, 0.93540077E+02,  0.31723213E+02, 0.89292258E+01, 0.21196877E+01,  0.43013488E+00, 0.75474799E-01, 0.11565621E-01,  0.15612002E-02, 0.18705605E-03, 0.20027866E-04,  0.19277885E-05, 0.16772276E-06];

  S2=[0.12000000E+02, 0.64800000E+02, 0.18514286E+03,  0.32400000E+03, 0.38173091E+03, 0.32300308E+03, 0.20566727E+03, 0.10207750E+03, 0.40592223E+02,  0.13221467E+02, 0.35916023E+01, 0.82606852E+00,   0.16293265E+00, 0.27862515E-01, 0.41703893E-02,  0.55091790E-03, 0.64704940E-04, 0.68008195E-05, 0.64341868E-06, 0.55082916E-07];

  R=abs(Z);
  RR=12.0;
  if(R-RR)<0
    N=max(1,min(20.,19.572/log(RR/R)));
  else
    N=20;
    println(" INBESS: RESULTS MAY BE INACCURATE ")
  end
  ZZ=Z/complex(RR,0.);
  ZZ=ZZ*ZZ;
  FINK0=-log(Z/complex(2.,0.))*(complex(1.,0.)+SER(S2,ZZ,1,N));
  FINK0=Z*(complex(S0,0.)+SER(S1,ZZ,1,N)+FINK0);

  return FINK0
end

function SER(S,ZZ,N1,N2)

  #  THIS FUNCTION COMPUTES THE VALUE OF A SERIES OF N2-N1+1
  #  TERMS CONSISTING OF INCREASING POWERS OF ZZ TIMES THE
  #  COEFFICIENTS IN ARRAY S

  y=complex(0.,0.);
  ZZZ=complex(1.,0.);
  i::Int64 = 1;
  for i=N1:N2
    ZZZ=ZZZ*ZZ;
    y=y+complex(S[i],0.)*ZZZ;
  end

  return y
end

function calc_solfundpot(x,y,xd,yd,nx,ny,tmp,k)
#Calcula as soluções fundamentais

r=sqrt((x-xd)^2+(y-yd)^2); # Raio (distância entre ponto fonte e
                           # ponto campo)
rx=(x-xd); # Componente x do raio
ry=(y-yd); # Componente y do raio
Tast=-1/(2*pi*k)*log(r); # Solução fundamental da temperatura
qast=1/(2*pi)*(rx*nx+ry*ny)/r^2; # Solução fundamental do fluxo
return Tast, qast
end



function  calc_solfund(x,y,xd,yd,nx,ny,CW,FR)
  #Calcula as soluções fundamentais
  r=sqrt((x-xd)^2+(y-yd)^2); # Raio (distância entre ponto fonte e
  # ponto campo)
  rx=(x-xd)/r; # Componente x do raio
  ry=(y-yd)/r; # Componente y do raio
  drdn=rx*nx+ry*ny;   #Componente do raio na direção normal
  ZR=real(FR*r/CW);
  Z=complex.(0.,ZR);
  F0C=SpecialFunctions.besselk(0,Z);
  F1C=SpecialFunctions.besselk(1,Z);

  qast=-(Z/r*drdn*F1C)/(2*pi); #Solução Fundamental da pressão acústica
  Tast=F0C/(2*pi);    #Solução Fundamental do fluxo de pressão acústica
  return Tast,qast
end

function trocacol(G,H,CDC)
  # Aplica as condições de contorno trocando as colunas das matrizes H e G
  ncdc = length(CDC); # número de linhas da matriz CDC
  A=H*1.0;
  B=G*1.0;
  for i=1:ncdc # Laço sobre as condições de contorno
    tipoCDC = CDC[i]; # Tipo da condição de contorno
    if tipoCDC == 0 # A temperatura é conhecida
      colunaA=-A[:,i]; # Coluna da matriz H que será trocada
      A[:,i]=-B[:,i]; # A matriz H recebe a coluna da matriz G
      B[:,i]=colunaA; # A mstriz G recebe a coluna da matriz H
    end
  end
  #  de contorno
  # b=B*valoresconhecidos; # vetor b
  return A, B
end
