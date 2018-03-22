function calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,CW,FR,fc,finc,qsi,w)
  # Calcula o potencial de velocidade nos pontos internos
  n_pint=length(PONTOS_int[:,1]); # Numero de pontos internos
  n_elem=length(phi); # Numero de elementos
  G_int = complex(zeros(n_pint,n_elem))
  H_int = complex(zeros(n_pint,n_elem))
  q = complex(zeros(n_pint,1))
  inc = complex(zeros(n_pint,1))

  for i=1:n_pint # La�o sobre os pontos internos
    x_fonte=PONTOS_int[i,2]; # Coordenada x do ponto fonte
    y_fonte=PONTOS_int[i,3]; # Coordenada y do ponto fonte
    #z_fonte=PONTOS_int[i,4]; # Coordenada z do ponto fonte
    for j=1:n_elem  #La�o sobre os elementos
      no1::Int64=ELEM[j,2]; # Primeiro n� geom�trico do elemento
      no2::Int64=ELEM[j,3]; # Segundo n� geom�trico do elemento

      x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
      y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1

      x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
      y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2

      G_int[i,j],H_int[i,j] = calcula_HeGns(x1,y1,x2,y2,x_fonte,y_fonte,CW,qsi,w,FR) # quando o ponto fonte nao pertence ao elemento
    end
        if fc[1,1] > 0
          q[i,1] = calc_q(x_fonte,y_fonte,fc,FR,CW)
          inc[i,1] = calc_inc(x_fonte,y_fonte,finc,FR,CW)
        else
          q[i,1], inc[i,1] =[0 0]
        end
  end

#  phi_pint = - (H_int*phi - G_int*qphi)
  phi_pint=-(H_int*phi-G_int*qphi-q-inc); # Vetor que contem o potencial de velocidade nos
  #      pontos internos
  return phi_pint
end

function calc_phi_pint_quad(PONTOS_int,NOS_GEO,ELEM,phi,qphi,CW,FR,fc,finc,qsi,w)
  # Calcula a temperatura nos pontos internos
  n_pint=length(PONTOS_int[:,1]); # Numero de pontos internos
  n_elem=length(phi); # Numero de elementos
  G_int = complex(zeros(n_pint,n_elem))
  H_int = complex(zeros(n_pint,n_elem))
  q = complex(zeros(n_pint,1))
  inc = complex(zeros(n_pint,1))

  for i=1:n_pint # La�o sobre os pontos internos
    x_fonte=PONTOS_int[i,2]; # Coordenada x do ponto fonte
    y_fonte=PONTOS_int[i,3]; # Coordenada y do ponto fonte
    #z_fonte=PONTOS_int[i,4]; # Coordenada z do ponto fonte
    for j=1:n_elem  #La�o sobre os elementos
      no1::Int64=ELEM[j,2]; # Primeiro n� geom�trico do elemento
      no2::Int64=ELEM[j,3]; # Segundo n� geom�trico do elemento
      no3::Int64=ELEM[j,4]; # Segundo n� geom�trico do elemento

      x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
      y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1

      x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
      y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2

      x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 2
      y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 2


      G_int[i,j],H_int[i,j] = calcula_HeGns_quad(x1,y1,x2,y2,x3,y3,x_fonte,y_fonte,CW,qsi,w,FR) # quando o ponto fonte nao pertence ao elemento
    end
        if fc[1,1] > 0
          q[i,1] = calc_q(x_fonte,y_fonte,fc,FR,CW)
          inc[i,1] = calc_inc(x_fonte,y_fonte,finc,FR,CW)
        else
          q[i,1], inc[i,1] =[0 0]
        end
  end

#  phi_pint = - (H_int*phi - G_int*qphi)
  phi_pint=-(H_int*phi-G_int*qphi-q-inc); # Vetor que contem a temperatura nos
  #      pontos internos
  return phi_pint
end
