function monta_GeH_linear(ELEM,NOS,k)
  #Geração dos pontos e pesos de Gauss
  npg=6;
  qsi,w = Gauss_Legendre(-1,1,npg)
  n_el = length(ELEM[:,1]);	# Numero total de elementos
  n_nos = length(NOS[:,1]);	# Numero total de nos
    # Inicializacao das matrizes H e G
  H = complex(zeros(n_nos,n_nos))
  G = complex(zeros(n_nos,2*n_el))
  for i = 1 : n_nos	# Percorre os pontos fontes
    # Coordenadas (xd,yd) dos pontos fontes
    xd = NOS[i,2]
    yd = NOS[i,3]
    for j = 1 : n_el	# Percorre os elementos
      # Numeracao dos dois nos do elemento j
      no1::Int64 = ELEM[j,2]
      no2::Int64 = ELEM[j,3]
      # Coordenadas dos dois nos (x1,y1,x2,y2)
      x1 = NOS[no1,2];	y1 = NOS[no1,3];
      x2 = NOS[no2,2];	y2 = NOS[no2,3];
      # Calculo das submatrizes h e g
      # Integracao nao singular
      # g,h = calc_gh_nsing(x1,y1,x2,y2,xd,yd,k);
      g,h = calcula_GeHns_linear(x1,y1,x2,y2,xd,yd,1,qsi,w,1)
      # println("diferença entre g e gii = ",abs(g-gii))

      if ((i == no1) || (i == no2))  # O no j pertence ao elemento j
        # gii=calc_g_sing(x1,y1,x2,y2,k); # Integracao singular com ponto fonte na extremidade do elemento qsi0=-1 ou qsi0=1)
      if(i==no1)
        qsitelles,Jtelles = telles(qsi,-1); # Calculando os pontos de Telles e o Jacobiano para o primeiro nó (qsi = -1)
        g,h = calcula_GeHns_linear(x1,y1,x2,y2,xd,yd,1,qsi,w,1)
        gii,h = calcula_GeHns_linear(x1,y1,x2,y2,xd,yd,1,qsitelles,w.*Jtelles,1)
        println("diferença entre g e gtelles = ",abs(g-gii))
        g=gii;
      else
        qsitelles,Jtelles = telles(qsi,1); # Calculando os pontos de Telles e o Jacobiano para o primeiro nó (qsi = 1)
        g,h = calcula_GeHns_linear(x1,y1,x2,y2,xd,yd,1,qsi,w,1)
        gii,h = calcula_GeHns_linear(x1,y1,x2,y2,xd,yd,1,qsitelles,w.*Jtelles,1)
        println("diferença entre g e gtelles = ",abs(g-gii))
        g=gii;
      end
      end
      for nolocal = 1 : 2
        noglobal::Int64 = ELEM[j,nolocal+1] #Indice da matriz global H
        H[i,noglobal] = H[i,noglobal] + h[nolocal]
      end
      G[i,2*j-1:2*j] = g;
    end
  end
  # Calculo dos termos da diagonal da matriz H (considera�ao de corpo a
  # temperatura constante)
  for m = 1 : n_nos
    H[m,m] = 0; #zera a diagonal principal
    for n = 1 : n_nos
      if n != m
        H[m,m] = H[m,m] - H[m,n]
      end
    end
  end
  return G,H
end

function calcula_GeHns_linear(x1,y1,x2,y2,xd,yd,CW,qsi,w,FR)
  #integração não singular
  n_pint=length(qsi); # Número de pontos de integração.
  g=complex(zeros(1,2)); # Inicializa o somatorio de g
  h=complex(zeros(1,2)); # Inicializa o somatorio de h
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

    k = CW/FR
    R=sqrt((x-xd)^2+(y-yd)^2); # Raio (dist�ncia entre ponto fonte e
    # ponto campo)
    Rx=(x-xd); # Componente x do raio
    Ry=(y-yd); # Componente y do raio
    Tast=-1/(2*pi*k)*log(R); # Solu��o fundamental da temperatura
    qast=1/(2*pi)*(Rx*nx+Ry*ny)/R^2; # Solu��o fundamental do fluxo
    # Tast,qast =calc_solfund(x,y,xd,yd,nx,ny,CW,FR); # Obtemos as soluções fundamentais

    h=h+qast*[N1 N2]*w[kk]*dgamadqsi; # Integral da matriz H
    g=g+Tast*[N1 N2]*w[kk]*dgamadqsi; # Integral da matriz G
  end
  return g,h
end

function calc_gh_nsing(x1,y1,x2,y2,xd,yd,k)

  # Pontos de Gauss
  qsi = [-0.86113631 -0.33998104 0.33998104 0.86113631];
  w = [0.34785485 0.65214515 0.65214515 0.34785485];

  npg=length(qsi); # N�mero de pontos de Gauss

  h=[0 0]; # Inicializa a matriz h do elemento
  g=[0 0]; # Inicializa a matriz g do elemento

  L=sqrt((x2-x1)^2+(y2-y1)^2); # Comprimento do elemento
  dgamadqsi=L/2; # Jacobiano

  sx=(x2-x1)/L; # Componente x do vetor tangente
  sy=(y2-y1)/L; # Componente y do vetor tangente
  nx=sy; # Componente x do vetor normal
  ny=-sx; # Componente y do vetor normal

  for kk=1:npg
    N1,N2=calc_fforma(qsi[kk]); # Calcula as fun��es de forma
    x=N1*x1+N2*x2; # Calcula a coordenada x do ponto de integra��o
    y=N1*y1+N2*y2; # Calcula a coordenada y do ponto de integra��o

    Tast,qast=calc_solfund(x,y,xd,yd,nx,ny,k,k); # Calcula as solu��es
    #  fundamentais
    h=h+qast*[N1 N2]*dgamadqsi*w[kk] # Integra��o da matriz h
    g=g+Tast*[N1 N2]*dgamadqsi*w[kk] # Integra��o da matriz g
  end
  return g,h
end

function aplica_CDC(G,H,CDC,ELEM)
  # Aplica as condi��es de contorno (troca as colunas das matrizes G e H)

  n_el=size(CDC,1);
  todos_valores=zeros(1,2*n_el);
  cont=0;

  # Cria a vari�vel T_PR que cont�m os n�s onde a temperatura � prescrita
  # (conhecida).
  # T_PR tem 5 colunas e o n�mero de linhas � igual ao n�mero de n�s onde a
  # temperatura � conhecida.
  # T_PR=[a1,a2,a3,a4,a5]
  # a1 = n�mero do elemento com temperatura prescrita ao qual este n�
  # pertence.
  # a2 = n�mero do n� com temperatura prescrita.
  # a3 = n�mero local do n� neste elemento.
  # a4: caso a temperatura tamb�m seja prescrita no segundo elemento a que
  # este n� pertence, ent�o, a4 conter� o n�mero deste elemento, caso
  # contr�rio, conter� zero.
  # a5: caso a4 seja diferente de zero, a5 conter� o n�mero local do n� no
  # segundo elemento, caso contr�rio, conter� zero.
  l,c = size(CDC)
  T_PR = zeros(l-1,c)
  H = zeros(l,l)
  G = zeros(l,l)
  troca = zeros(l,l)

  for el=1:n_el # for sobre os elementos para criar T_PR
    for no=1:2 # for sobre os n�s do elemento el
      no_global=ELEM[el,no+1]; # n�mero global do n�
      tipoCDC=CDC[el,2*no];  # tipo da condi��o de contorno
      valorCDC=CDC[el,2*no+1]; # valor da condi��o de contorno
      todos_valores[2*el-2+no]=valorCDC; # armazerna o valor da condi��o
      # de contorno no vetor todos_valores
      if(tipoCDC==0) # se a temperatura � conhecida
        compartilha=0; # por enquanto n�o se sabe se a tempertura �
        # tamb�m conhecida no segundo elemento a que este n� pertence
        i=1;
        while (!Bool(compartilha)&&i<=cont&&cont>0) # verifica se o n� global
          # j� est� presente em T_PR. Quando ele encontra,
          # compartilha se torna igual a 1 e o while p�ra.
          if(no_global==T_PR[i,2]) # Se sim, o n� global j� est�
            # presente em T_PR. As colunas 4 e 5 s�o preenchidas e
            # compartilha se torna igual a 1.
            compartilha=1;
            T_PR[i,4]=el;
            T_PR[i,5]=no;
          end
          i=i+1;
        end
        if(!Bool(compartilha)) # Se compartilha continua zero, ent�o o n�
          # ainda n�o foi inserido em T_PR. Neste caso, as tr�s
          # primeiras colunas de T_PR s�o preenchidas.
          cont=cont+1;
          T_PR[cont,1]=el;
          T_PR[cont,2]=no_global;
          T_PR[cont,3]=no;
        end
      end
    end
  end

  n_temp_pr = length(T_PR[:,1]); # N�mero de n�s com temperatura conhecidas
  for i=1 : n_temp_pr # for sobre os n�s com temperatura conhecida
    i_el=T_PR[i,1]; # primeiro elemento que cont�m o n�
    i_no=T_PR[i,2]; # n�mero do n� com temperatura conhecida
    i_no_loc=T_PR[i,3]; # n�mero local do n� neste elemento
    ind_H::Int64=i_no; # �ndice da coluna da matriz H que ser� trocada
    ind_G::Int64 =2*i_el+i_no_loc-2; # �ndice da coluna da matriz G que ser�
    #  trocada
    troca = G[:,ind_G]; # armazena a coluna de G que ser� trocada
    G[:,ind_G] = -H[:,ind_H]; # substitui a coluna de H na de G
    H[:,ind_H] = -troca; # Substitui a de G na de H
    if(T_PR[i,4]!=0) # Se a temperatura tamb�m � conhecida no segundo
      #   elemento
      i_el=T_PR[i,4]; # N�mero do segundo elemento
      i_no_loc=T_PR[i,5]; # n�mero local deste n� no segundo elemento
      ind_G=2*i_el+i_no_loc-2; # �ndice da coluna G que ser� atribu�do
      # zero
      H[:,ind_H]=H[:,ind_H]-G[:,ind_G]; # soma o valor da coluna G na
      # coluna H
      G[:,ind_G]=0; # atribui zero na coluna G
    end
  end
  b=G*todos_valores'; # C�lculo do vetor b
  A=H;
  return A,b,T_PR
end
