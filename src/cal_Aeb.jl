function cal_Aeb(b1,b2,arg)
  NOS,NOS_GEO,ELEM,CDC,CW,FR,qsi,w=arg
  nelem::Int64=size(ELEM)[1]; # Numero de elementos (numero de linhas da
  #  matriz ELEM)
  qsitelles,Jtelles = telles(qsi,0); # Calculando os pontos de Telles e o Jacobiano
  G=complex(zeros(length(b1),length(b2)));
  H=complex(zeros(length(b1),length(b2)));
  q=zeros(length(b1),1);
  ci=0
  for i in b1 # Laco sobre os pontos fontes
    ci+=1
    xd=NOS[i,2]; # Coordenada x do ponto fonte
    yd=NOS[i,3]; # Coordenada y do ponto fonte
    cj=0
    for j in b2 # Laco sobre os elementos
      cj+=1
      noi::Int64=ELEM[j,2]; # Ponto inicial do elemento
      nof::Int64=ELEM[j,3]; # Ponto final do elemento
      x1=NOS_GEO[noi,2]; # Coordenada x do ponto inicial do elemento
      x2=NOS_GEO[nof,2]; # Coordenada x do ponto final do elemento
      y1=NOS_GEO[noi,3]; # Coordenada y do ponto inicial do elemento
      y2=NOS_GEO[nof,3];  # Coordenada y do ponto final do elemento
      if i==j # O ponto fonte pertence ao elemento
        g,h = calcula_HeGs(x1,y1,x2,y2,CW,FR);
        gtelles,htelles = calcula_HeGns(x1,y1,x2,y2,xd,yd,CW,qsitelles,w.*Jtelles,FR);
        htelles = htelles + 0.5
        # println("Diferença entre g e gtelles = ", abs(g-gtelles))
        # println("Diferença entre h e htelles = ", abs(h-htelles))
      else # O ponto fonte n�o pertence ao elemento
        g,h = calcula_HeGns(x1,y1,x2,y2,xd,yd,CW,qsitelles,w.*Jtelles,FR);
      end
      if CDC[j,2]==0
        G[ci,cj] = -h
        H[ci,cj] = -g
      else
        G[ci,cj] = g
        H[ci,cj] = h
      end
    end
  end
  return H,G
end

function cal_Aeb_linear(b1,b2,arg)
  NOS,NOS_GEO,ELEM,CDC,CW,FR,qsi,w=arg
  nelem::Int64=size(ELEM)[1]; # Numero de elementos (numero de linhas da
  #  matriz ELEM)
  G=complex(zeros(length(b1),length(b2)));
  H=complex(zeros(length(b1),length(b2)));
  q=zeros(length(b1),1);
  ci=0
  for i in b1 # Laco sobre os pontos fontes
    ci+=1
    xd=NOS[i,2]; # Coordenada x do ponto fonte
    yd=NOS[i,3]; # Coordenada y do ponto fonte
    cj=0
    for j in b2 # Laco sobre os elementos
      cj+=1
      noi::Int64=ELEM[j,2]; # Ponto inicial do elemento
      nof::Int64=ELEM[j,3]; # Ponto final do elemento
      x1=NOS_GEO[noi,2]; # Coordenada x do ponto inicial do elemento
      x2=NOS_GEO[nof,2]; # Coordenada x do ponto final do elemento
      y1=NOS_GEO[noi,3]; # Coordenada y do ponto inicial do elemento
      y2=NOS_GEO[nof,3];  # Coordenada y do ponto final do elemento
      if i==no1 # O ponto fonte pertence ao elemento
        # g,h = calcula_HeGns_quad(x1,y1,x2,y2,x3,y3,CW,FR);
        qsitelles,Jtelles = telles(qsi,-1); # Calculando os pontos de Telles e o Jacobiano
        g,h = calcula_HeGns_linear(x1,y1,x2,y2,xd,yd,CW,qsitelles,w.*Jtelles,FR)
      elseif i==no2 # O ponto fonte pertence ao elemento
        qsitelles,Jtelles = telles(qsi,1); # Calculando os pontos de Telles e o Jacobiano
        g,h = calcula_HeGns_linear(x1,y1,x2,y2,xd,yd,CW,qsitelles,w.*Jtelles,FR)
      end # O ponto fonte não pertence ao elemento

      g,h = calcula_HeGns_linear(x1,y1,x2,y2,xd,yd,CW,qsi,w,FR)

      # Determinar se o nó pertence ao elemento
      for nolocal = 1 : 2
        noglobal = ELEM[j,nolocal+1]; #�ndice da matriz global H
        H[i,noglobal] = H[i,noglobal] + h[nolocal];
      end
      G[i,2*j-1:2*j] = g

      if CDC[j,2]==0
        G[ci,cj] = -h
        H[ci,cj] = -g
      else
        G[ci,cj] = g
        H[ci,cj] = h
      end
    end
  end
  return H,G
end

function calcula_HeGns_linear(x1,y1,x2,y2,xd,yd,CW,qsi,w,FR)
  #integração não singular
  n_pint=length(qsi); # Número de pontos de integração.
  g=complex(zeros(2)); # Inicializa o somatorio de g
  h=complex(zeros(2)); # Inicializa o somatorio de h
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
    h=h+[N1 N2].*qast*w[kk]*dgamadqsi; # Integral da matriz H
    g=g+[N1 N2].*Tast*w[kk]*dgamadqsi; # Integral da matriz G
  end
  return g,h
end


function cal_Aeb_quad(b1,b2,arg)
  NOS,NOS_GEO,ELEM,CDC,CW,FR,qsi,w=arg
  nelem::Int64=size(ELEM)[1]; # Numero de elementos (numero de linhas da
  #  matriz ELEM)
  G=complex(zeros(length(b1),length(b1)));
  H=complex(zeros(length(b1),length(b1)));
  ci=0
  for i in b1 # Laco sobre os pontos fontes
    ci+=1
    xd=NOS[i,2]; # Coordenada x do ponto fonte
    yd=NOS[i,3]; # Coordenada y do ponto fonte
    cj=0
    for j in b2 # Laco sobre os elementos
      cj+=1
      no1::Int64=ELEM[j,2]; # Ponto inicial do elemento
      no2::Int64=ELEM[j,3]; # Ponto intermediário do elemento
      no3::Int64=ELEM[j,4]; # Ponto final do elemento
      x1=NOS_GEO[no1,2]; # Coordenada x do ponto inicial do elemento
      x2=NOS_GEO[no2,2]; # Coordenada x do ponto intermediário do elemento
      x3=NOS_GEO[no3,2]; # Coordenada x do ponto final do elemento
      y1=NOS_GEO[no1,3]; # Coordenada x do ponto inicial do elemento
      y2=NOS_GEO[no2,3]; # Coordenada x do ponto intermediário do elemento
      y3=NOS_GEO[no3,3]; # Coordenada x do ponto final do elemento
      if i==no1 # O ponto fonte pertence ao elemento
        # g,h = calcula_HeGns_quad(x1,y1,x2,y2,x3,y3,CW,FR);
        qsitelles,Jtelles = telles(qsi,-1); # Calculando os pontos de Telles e o Jacobiano
        g,h = calcula_HeGns_quad(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsitelles,w.*Jtelles,FR);
      elseif i==no2 # O ponto fonte pertence ao elemento
        qsitelles,Jtelles = telles(qsi,0); # Calculando os pontos de Telles e o Jacobiano
        g,h = calcula_HeGns_quad(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsitelles,w.*Jtelles,FR);
      elseif i==no3
        qsitelles,Jtelles = telles(qsi,1); # Calculando os pontos de Telles e o Jacobiano
        g,h = calcula_HeGns_quad(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsitelles,w.*Jtelles,FR);
      end # O ponto fonte não pertence ao elemento
      g,h = calcula_HeGns_quad(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsi,w,FR);
      for nolocal = 1 : 3
          noglobal::Int64 = ELEM[j,nolocal+1]; # Indice da matriz global H
          H[i,noglobal] = H[i,noglobal] + h[nolocal];
        #   G[i,noglobal] = g[nolocal]
      end
      G[i,3*j-2:3*j] = g

      # if CDC[j,2]==0
      #   G[ci,cj] = -h
      #   H[ci,cj] = -g
      # else
      #   G[ci,cj] = g
      #   H[ci,cj] = h
      # end
    end
  end
  return H,G
end

function calcula_HeGns_quad(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsi,w,FR);
  #integração não singular
  n_pint=length(qsi); # Número de pontos de integração.
  g=complex(zeros(3)); # Inicializa o somatorio de g
  h=complex(zeros(3)); # Inicializa o somatorio de h
  N = zeros(n_pint,3)
  dN = zeros(n_pint,3)
  dgamadqsi = zeros(n_pint,1)
  phiast =complex(zeros(n_pint,1))
  qast = complex(zeros(n_pint,1))
  for kk=1:n_pint # Laço sobre os pontos de integração
    N[kk,:] =calc_fforma_quad(qsi[kk]); # Calcula as funções de forma
    dN[kk,:] = calc_dfforma_quad(qsi[kk])
    xx=N[kk,1]*x1+N[kk,2]*x2+N[kk,3]*x3; # Calcula a coordenada x do ponto de integração
    yy=N[kk,1]*y1+N[kk,2]*y2+N[kk,3]*y3; # Calcula a coordenada y do ponto de integração
    dx=dN[kk,1]*x1+dN[kk,2]*x2+dN[kk,3]*x3; # Calcula a coordenada x do ponto de integração
    dy=dN[kk,1]*y1+dN[kk,2]*y2+dN[kk,3]*y3; # Calcula a coordenada y do ponto de integração
    dgamadqsi[kk]=cal_Jacobiano([x1 x2 x3],[y1 y2 y3],dN[kk,:]); # Jacobiano
    sx=dx/dgamadqsi[kk]; # Componente x do vetor tangente
    sy=dy/dgamadqsi[kk]; # Componente y do vetor tangente
    nx=sy; # Componente x do vetor normal
    ny=-sx; # Componente y do vetor normal
    phiast[kk],qast[kk] =calc_solfund(xx,yy,xd,yd,nx,ny,CW,FR); # Obtemos as soluções fundamentais
    #phiast[kk],qast[kk] = [1 1]
    h=h+N[kk,:]*phiast[kk]*w[kk]*dgamadqsi[kk]; # Integral da matriz H
    g=g+N[kk,:]*qast[kk]*w[kk]*dgamadqsi[kk]; # Integral da matriz G
  end
  return g,h
end

# Só lembrando que podemos utilizar a cal_Aeb_quad
# e calcula_HeGns_quad apenas modificando a função de forma.
function cal_Aeb_Bezier(b1,b2,arg)
  NOS,NOS_GEO,ELEM,CDC,CW,FR,qsi,w=arg
  nelem::Int64=size(ELEM)[1]; # Numero de elementos (numero de linhas da
  #  matriz ELEM)
  G=complex(zeros(length(b1),length(b2)));
  H=complex(zeros(length(b1),length(b2)));
  ci=0
  for i in b1 # Laco sobre os pontos fontes
    ci+=1
    xd=NOS[i,2]; # Coordenada x do ponto fonte
    yd=NOS[i,3]; # Coordenada y do ponto fonte
    cj=0
    for j in b2 # Laco sobre os elementos
      cj+=1
      no1::Int64=ELEM[j,2]; # Ponto inicial do elemento
      no2::Int64=ELEM[j,3]; # Ponto intermediário do elemento
      no3::Int64=ELEM[j,4]; # Ponto final do elemento
      x1=NOS_GEO[no1,2]; # Coordenada x do ponto inicial do elemento
      x2=NOS_GEO[no2,2]; # Coordenada x do ponto intermediário do elemento
      x3=NOS_GEO[no3,2]; # Coordenada x do ponto final do elemento
      y1=NOS_GEO[no1,3]; # Coordenada x do ponto inicial do elemento
      y2=NOS_GEO[no2,3]; # Coordenada x do ponto intermediário do elemento
      y3=NOS_GEO[no3,3]; # Coordenada x do ponto final do elemento
      if i==j # O ponto fonte pertence ao elemento
        # g,h = calcula_HeGns_quad(x1,y1,x2,y2,x3,y3,CW,FR);
        g,h = calcula_HeGns_quad(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsi,w,FR);
      else # O ponto fonte n�o pertence ao elemento
        g,h = calcula_HeGns_quad(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsi,w,FR);
      end
      if CDC[j,2]==0
        G[ci,cj] = -h
        H[ci,cj] = -g
      else
        G[ci,cj] = g
        H[ci,cj] = h
      end
    end
  end
  return H,G
end

function calcula_HeGns_Bezier(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsi,w,FR);
  #integração não singular
  n_pint=length(qsi); # Número de pontos de integração.
  g=complex(0,0); # Inicializa o somatorio de g
  h=complex(0,0); # Inicializa o somatorio de h
  N = zeros(size(qsi,2),3)
  dN = zeros(size(qsi,2),3)
  dgamadqsi = zeros(size(qsi,2),1)
  phiast =complex(zeros(size(qsi,2),1))
  qast = complex(zeros(size(qsi,2),1))
  for kk=1:n_pint # Laço sobre os pontos de integração
    N[kk,:] =calc_fforma_Bezier(qsi[kk]); # Calcula as funções de forma
    dN[kk,:] = calc_dfforma_Bezier(qsi[kk])
    xx=N[kk,1]*x1+N[kk,2]*x2+N[kk,3]*x3; # Calcula a coordenada x do ponto de integração
    yy=N[kk,1]*y1+N[kk,2]*y2+N[kk,3]*y3; # Calcula a coordenada y do ponto de integração
    dx=dN[kk,1]*x1+dN[kk,2]*x2+dN[kk,3]*x3; # Calcula a coordenada x do ponto de integração
    dy=dN[kk,1]*y1+dN[kk,2]*y2+dN[kk,3]*y3; # Calcula a coordenada y do ponto de integração
    dgamadqsi[kk]=cal_Jacobiano([x1 x2 x3],[y1 y2 y3],dN[kk,:]); # Jacobiano
    sx=dx/dgamadqsi[kk]; # Componente x do vetor tangente
    sy=dy/dgamadqsi[kk]; # Componente y do vetor tangente
    nx=sy; # Componente x do vetor normal
    ny=-sx; # Componente y do vetor normal
    phiast[kk],qast[kk] =calc_solfund(xx,yy,xd,yd,nx,ny,CW,FR); # Obtemos as soluções fundamentais
    #phiast[kk],qast[kk] = [1 1]
    h=h+phiast[kk]*w[kk]*dgamadqsi[kk]; # Integral da matriz H
    g=g+qast[kk]*w[kk]*dgamadqsi[kk]; # Integral da matriz G
  end
  return g,h
end
