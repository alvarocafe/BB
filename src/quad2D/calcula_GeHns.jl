function calcula_GeHns(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsi,w,FR)
# Non singular integration
  n_pint=size(qsi,1); # Number of integration points
  h = complex(zeros(3)) # Initializes the matrix H terms
  g = complex(zeros(3)) # Initializes the matrix G terms

  for kk=1:n_pint # Loop over the integration points
  N=calc_fforma(qsi[kk]); # Evaluates the shape functions
  dN= calc_dfforma(qsi[kk]) # Evaluates the derivative of the shape functions
  dgamadqsi=cal_Jacobiano([x1 x2 x3],[y1 y2 y3],dN); # Evaluates the Jacobian
  xx=N[1]*x1+N[2]*x2+N[3]*x3; # x coordinate of the integration point
  yy=N[1]*y1+N[2]*y2+N[3]*y3; # y coordinate of the integration point
  dx=dN[1]*x1+dN[2]*x2+dN[3]*x3; # x coordinate derivative of the curve at the integration point
  dy=dN[1]*y1+dN[2]*y2+dN[3]*y3; # x coordinate derivative of the curve at the integration point
  sx=dx/dgamadqsi # x component of the tangent vector
  sy=dy/dgamadqsi # y component of the tangent vector
  nx=sy; # x component of the normal vector
  ny=-sx; # y component of the normal vector
#  phiast,qast =calc_solfund(xx,yy,xd,yd,nx,ny,CW,FR); # Evaluation of the fundamental solutions
	phiast = 1	#uncomment these lines to obtain the length of the boundary
	qast = 1
  g=g+N'*phiast*w[kk]*dgamadqsi # Evaluation of integral for g term
  h=h+N'*qast*w[kk]*dgamadqsi # Evaluation of integral for h term
  end
  return g,h
end

