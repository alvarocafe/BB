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
