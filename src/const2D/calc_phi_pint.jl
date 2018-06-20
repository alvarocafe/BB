function calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,finc,qsi,w,k)
# Evaluates the velocity potential at internal (or external) points
  n_pint=length(PONTOS_int[:,1]); # Number of internal points
  n_elem=length(phi); # Number of elements
  G_int = complex(zeros(n_pint,n_elem))	# Allocates the matrix for the influence of the flux at the internal points
  H_int = complex(zeros(n_pint,n_elem)) # Allocates the matrix for the influence of the velocity potential at the internal points
  q = complex(zeros(n_pint,1))	# Allocates the array for the influence of concentrated acoustic sources
  inc = complex(zeros(n_pint,1))	# Allocates the array for the influence of incident plane waves

  for i=1:n_pint # Loop over the internal points
    x_fonte=PONTOS_int[i,2]; # x coordinate of the internal point
    y_fonte=PONTOS_int[i,3]; # y coordinate of the internal point
    for j=1:n_elem  # Loop over the elements
      no1::Int64=ELEM[j,2]; # First point of the element
      no2::Int64=ELEM[j,3]; # Second point of the element

      x1=NOS_GEO[no1,2]; # x coordinate of the first point
      y1=NOS_GEO[no1,3]; # y coordinate of the first point

      x2=NOS_GEO[no2,2]; # x coordinate of the second point
      y2=NOS_GEO[no2,3]; # y coordinate of the second point

      G_int[i,j],H_int[i,j] = calcula_GeHns(x1,y1,x2,y2,x_fonte,y_fonte,qsi,w,k) # Non singular integration
    end
        if fc[1,1] > 0 	# If the flag is greater than zero, there are acoustic sources and/or incident plane waves
          q[i,1] = calc_q(x_fonte,y_fonte,fc,k) # Evaluates the influence from acoustic concentrated sources
          inc[i,1] = calc_inc(x_fonte,y_fonte,finc,k) # Evaluates the influence from acoustic incident plane waves
        else
          q[i,1], inc[i,1] =[0 0]	# There is no influence
        end
  end
#  phi_pint=-(H_int*phi-G_int*qphi-q-inc); # Evaluates the velocity potential for the internal points
  phi_pint=-(H_int*phi-G_int*qphi); # Evaluates the velocity potential for the internal points
 
  return phi_pint
end

function calc_phi_pintpot(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,qsi,w,k)
# Evaluates the velocity potential at internal (or external) points
  n_pint=length(PONTOS_int[:,1]); # Number of internal points
  n_elem=length(phi); # Number of elements
  G_int = complex(zeros(n_pint,n_elem))	# Allocates the matrix for the influence of the flux at the internal points
  H_int = complex(zeros(n_pint,n_elem)) # Allocates the matrix for the influence of the velocity potential at the internal points
  q = complex(zeros(n_pint,1))	# Allocates the array for the influence of concentrated acoustic sources
  inc = complex(zeros(n_pint,1))	# Allocates the array for the influence of incident plane waves

  for i=1:n_pint # Loop over the internal points
    x_fonte=PONTOS_int[i,2]; # x coordinate of the internal point
    y_fonte=PONTOS_int[i,3]; # y coordinate of the internal point
    for j=1:n_elem  # Loop over the elements
      no1::Int64=ELEM[j,2]; # First point of the element
      no2::Int64=ELEM[j,3]; # Second point of the element

      x1=NOS_GEO[no1,2]; # x coordinate of the first point
      y1=NOS_GEO[no1,3]; # y coordinate of the first point

      x2=NOS_GEO[no2,2]; # x coordinate of the second point
      y2=NOS_GEO[no2,3]; # y coordinate of the second point

      G_int[i,j],H_int[i,j] = calcula_GeHnspot(x1,y1,x2,y2,x_fonte,y_fonte,qsi,w,k) # Non singular integration
    end
        if fc[1,1] > 0 	# If the flag is greater than zero, there are acoustic sources and/or incident plane waves
          q[i,1] = calc_q(x_fonte,y_fonte,fc,k) # Evaluates the influence from acoustic concentrated sources
        else
          q[i,1] =0	# There is no influence
        end
  end
#  phi_pint=-(H_int*phi-G_int*qphi-q-inc); # Evaluates the velocity potential for the internal points
  phi_pint=(H_int*phi-G_int*qphi); # Evaluates the velocity potential for the internal points
 
  return phi_pint
end

function domain_field(xfield,y,T,q,node,dnorm,kmat)
# -------------------------------------------------------------------------------
#  Fucntion used to evaluate the temperature/heat flux using the direct
#     evaluation of the integral representation
#
# -------------------------------------------------------------------------------
if(isempty(xfield))
  f=zeros(1)
  fx=zeros(1)
  fy=zeros(1)
else
    pi2 = π*2
    nfield=size(xfield,2)
    f=zeros(nfield)
    fx=zeros(nfield)
    fy=zeros(nfield)
    n=size(y,2)
    al= sqrt.((y[1,vec(node[2,:])]-y[1,vec(node[1,:])]).^2+ (y[2,vec(node[2,:])]-y[2,vec(node[1,:])]).^2)

    for j=1:n
        # Comprimento do elemento
        # Dist�ncias entre o ponto interno e os extremos dos elementos
        x11 = y[1,node[1,j]] - xfield[1,:]
        x21 = y[2,node[1,j]] - xfield[2,:]
        x12 = y[1,node[2,j]] - xfield[1,:]
        x22 = y[2,node[2,j]] - xfield[2,:]
        r1 =  sqrt.(x11.^2 + x21.^2); # Dist�ncia para o in�cio do elemento
        r2 =  sqrt.(x12.^2 + x22.^2); # Dist�ncia para o final do elemento

        # Proje��o do vetor dist�ncia no vetor normal ao elemento
        d  =  x11.*dnorm[1,j] + x21.*dnorm[2,j]; # Figura A.1, p�gina 178
        t1 = -x11.*dnorm[2,j] + x21.*dnorm[1,j]; # Dist�ncia T1 da figura A.1
        t2 = -x12.*dnorm[2,j] + x22.*dnorm[1,j]; # Dist�ncia T2 da figura A.1
        ds = abs.(d)
        dtheta = atan2.(ds.*al[j],ds.^2+t1.*t2)

        # Equa��o (A.5) do livro do Liu: elementos da matriz G
        g = -(-dtheta.*ds + al[j] + t1.*real(log.(r1))-t2.*real(log.(r2)))/(pi2*kmat)
        # Equa��o (A.7) com nx=1, ny=0, tx=0, ty=1.
        kkx = (dtheta.*dnorm[1,j] - real(log.(r2./r1)).*dnorm[2,j])/(pi2*kmat)
        # Equa��o(A.7) do livro do Liu com nx=0, ny=1, tx=-1, ty=0.
        kky = (dtheta.*dnorm[2,j] + real(log.(r2./r1)).*dnorm[1,j])/(pi2*kmat)
        hhx = -(-(t2./r2.^2-t1./r1.^2).*dnorm[1,j] - d.*(1./r2.^2-1./r1.^2).*dnorm[2,j])/pi2; # Equa��o (A.8) com nx=1
        # ny=0; tx=0; ty=1.
        hhy = -(-(t2./r2.^2-t1./r1.^2).*dnorm[2,j] + d.*(1./r2.^2-1./r1.^2).*dnorm[1,j])/pi2; # Equa��o (A.8) com nx=0
        #     if(d<=0)
        #         dtheta = -dtheta
        #     end
        dtheta = dtheta.*sign.(d)
        h = -dtheta/pi2; # Equa��o (A.6): Elementos da matriz 
        f = f + g*q[j] - h*T[j]; # Integral (2.12) com o termo de
        # dom�nio igual a zero.
        fx = fx + kkx*q[j] - hhx*T[j]; # Integral (2.12) com o termo de
        # dom�nio igual a zero.
        fy = fy + kky*q[j] - hhy*T[j]; # Integral (2.12) com o termo de
        # dom�nio igual a zero.
    end
end
return f,fx,fy
end
