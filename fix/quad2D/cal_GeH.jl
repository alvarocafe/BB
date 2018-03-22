function cal_GeH(NOS,ELEM,CW,FR,fc,qsi,w)
# Evaluates the G and H matrices for the linear system H phi = G q, where phi is a vector containing the values of the velocity potential and q is a vector containing the values of the flux at the boundary of the problem.
 
  nelem::Int8=size(ELEM)[1]; # Number of elements
  nnos::Int8=size(NOS)[1]; # Number of nodes
  G=complex(zeros(nnos,nnos)); 	# Allocates matrix G
  H=complex(zeros(nnos,nnos));	# Allocates matrix H
  q=zeros(nnos,1);  # Influence from concentrated sources
  inc = zeros(nnos,1);  # Influence from incident plane waves
  qsitelles1,Jtelles1 = telles(qsi,-1); # Evaluates the Telles' points and Jacobian if the field point corresponds to the first node of the element
  qsitelles2,Jtelles2 = telles(qsi,0); # Evaluates the Telles' points and Jacobian if the field point corresponds to the second node of the element
  qsitelles3,Jtelles3 = telles(qsi,1); # Evaluates the Telles' points and Jacobian if the field point corresponds to the third node of the element

  for i=1:nnos # Loop over the source points
    xd=NOS[i,2]; # x coordinate of the source point
    yd=NOS[i,3]; # y coordinate of the source point
    for j=1:nelem # Loop over the elements
      no1::Int8=ELEM[j,2]; # First point of the element
      no2::Int8=ELEM[j,3]; # Second point of the element
      no3::Int8=ELEM[j,4]; # Second point of the element
      x1=NOS[no1,2]; # x coordinate of the first point of the element
      x2=NOS[no2,2]; # x coordinate of the second point of the element
      x3=NOS[no3,2]; # x coordinate of the second point of the element
      y1=NOS[no1,3]; # y coordinate of the first point of the element
      y2=NOS[no2,3];  # y coordinate of the second point of the element
      y3=NOS[no3,3];  # y coordinate of the second point of the element
      if i==no1 # The source point belongs to the element
        g, h = calcula_GeHns(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsitelles1,w.*Jtelles1,FR);	# Singular integration using the Telles transformation
        #h = h - 0.5	# Adding the jump term
        elseif i==no2 # The source point doesnt belong to the element
	g, h = calcula_GeHns(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsitelles2,w.*Jtelles2,FR);	# Singular integration using the Telles transformation
        #h = h - 0.5	# Adding the jump term
	elseif i==no3
	g, h = calcula_GeHns(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsitelles3,w.*Jtelles3,FR);	# Singular integration using the Telles transformation      
        #h = h - 0.5	# Adding the jump term
	else
	g, h = calcula_GeHns(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsi,w,FR);	#Regular integration using Gaussian quadrature
	end
	println("Tipo e tamanho de g:\n tipo=",typeof(g),",\n tamanho=",size(g),",\n valor g=",g)
	println("Tipo e tamanho de h:\n tipo=",typeof(h),",\n tamanho=",size(h),",\n valor h=",h)
	# allocate evaluated integrals into the global H and G matrices
 	for nolocal=1:3
            noglobal::Int8=ELEM[j,nolocal+1]; #Global index of matrix H
            H[i,noglobal] = H[i,noglobal] + h[nolocal];
        end
	println("Elemento = $j")
            G[i,3*j-2:3*j] = g
    end
  end
  # Rigid body consideration for the evaluation of the singular terms of matrix H
for m = 1 : nnos
    H[m,m] = 0; # let the diagonal be zero
    for n = 1 : nnos	# Loop over the source points
        if n != m	# If it's not a diagonal term
            H[m,m] = H[m,m] - H[m,n]; # Subtracts the term from the diagonal
        end
    end
end

  return G,H
end
