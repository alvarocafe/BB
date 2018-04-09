function aplica_CDC(G,H,CDC)
  # Aplica as condições de contorno trocando as colunas das matrizes H e G
  ncdc = length(CDC[:,1]); # número de linhas da matriz CDC
  A=H*1.0;
  B=G*1.0;
  for i=1:ncdc # Laço sobre as condições de contorno
    tipoCDC = CDC[i,2]; # Tipo da condição de contorno
    if tipoCDC == 0 # A temperatura é conhecida
      colunaA=-A[:,i]; # Coluna da matriz H que será trocada
      A[:,i]=-B[:,i]; # A matriz H recebe a coluna da matriz G
      B[:,i]=colunaA; # A mstriz G recebe a coluna da matriz H
    end
  end
  valoresconhecidos=CDC[:,3]; # Valores das condições
  #  de contorno
  b=B*valoresconhecidos; # vetor b
  return A, b
end

function monta_Teq(CDC,x)
  # Separa fluxo e temperatura

  # ncdc = n�mero de linhas da matriz CDC
  # T = vetor que cont�m as temperaturas nos n�s
  # q = vetor que cont�m o fluxo nos n�s

  ncdc = length(CDC[:,1]);
  nnos = length(x)
  T = complex(zeros(nnos,1))
  q = complex(zeros(nnos,1))
  for i=1:ncdc # Laco sobre as condicoes de contorno
    tipoCDC=CDC[i,2]; # Tipo da condi��o de contorno
    valorCDC=CDC[i,3]; # Valor da condi��o de contorno
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
