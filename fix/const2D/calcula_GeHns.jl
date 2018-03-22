function calcula_GeHns(x1,y1,x2,y2,xd,yd,qsi,w,k)
# Non singular integration
  n_pint=length(qsi); # Number of integration points
  g=complex(0,0); # Start the sum for g
  h=complex(0,0); # Start the sum for h
  L=sqrt((x2-x1)^2+(y2-y1)^2); # Element length
  dgamadqsi=L/2; # Jacobian

  sx=(x2-x1)/L; # x component of the tangent vector
  sy=(y2-y1)/L; # y component of the tangent vector
  nx=sy; # x component of the normal vector
  ny=-sx; # y component of the normal vector

  for kk=1:n_pint # Loop over the integration points
    N1,N2 =calc_fforma(qsi[kk]); # Evaluates the shape functions
    x=N1*x1+N2*x2; # Evaluates the x coordinate of the integration point
    y=N1*y1+N2*y2; # Evaluates the y coordinate of the integration point

    Tast,qast =calc_solfund(x,y,xd,yd,nx,ny,k); # Evaluate the fundamental solutions
    #Tast = 1; qast = 1; #To test the integration, let the fundamental solutions be 1. The result of the sum of each line of matrices H and G should be the perimeter of the geometric model.
    h=h+qast*w[kk]*dgamadqsi; # Integration of matrix H
    g=g+Tast*w[kk]*dgamadqsi; # Integration of matrix G
  end
  return g,h
end

function calcula_GeHnspot(x1,y1,x2,y2,xd,yd,qsi,w,k)
#integração não singular
n_pint=length(qsi); # Número de pontos de integração (comprimento do
                  #    vetor qsi)
G=0; # Inicializa o somatorio de g
H=0; # Inicializa o somatorio de h
for kk=1:n_pint # Laço sobre os pontos de integração
    N1,N2=calc_fforma(qsi[kk]); # Calcula as funções de forma
    dN1dqsi= -0.5
    dN2dqsi=0.5 # Calcula as derivadas das
                                                     #    funções de forma
    x=N1*x1+N2*x2; # Calcula a coordenada x do ponto de integração
    y=N1*y1+N2*y2; # Calcula a coordenada y do ponto de integração

    dxdqsi=dN1dqsi*x1+dN2dqsi*x2;
    dydqsi=dN1dqsi*y1+dN2dqsi*y2;
    dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);

    sx=dxdqsi/dgamadqsi; # Componente x do vetor tangente
    sy=dydqsi/dgamadqsi; # Componente y do vetor tangente
    nx=sy; # Componente x do vetor normal
    ny=-sx; # Componente y do vetor normal

    Tast,qast=calc_solfundpot(x,y,xd,yd,nx,ny,k); # Calcula as soluções fundamentais
    H=H+qast*dgamadqsi*w[kk]; # Integral da matriz H
    G=G+Tast*dgamadqsi*w[kk]; # Integral da matriz G
end
return G, H
end

