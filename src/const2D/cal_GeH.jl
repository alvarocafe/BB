function cal_GeH(NOS,NOS_GEO,ELEM,k,fc,qsi,w)
# Evaluates the G and H matrices for the linear system H phi = G q, where phi is a vector containing the values of the velocity potential and q is a vector containing the values of the flux at the boundary of the problem.
 
  nelem::Int64=size(ELEM)[1]; # Number of elements
  nnos::Int64=nelem; # Number of nodes
  G=complex(zeros(nnos,nnos)); 	# Allocates matrix G
  H=complex(zeros(nnos,nnos));	# Allocates matrix H
  q=zeros(nnos,1);  # Influence from concentrated sources
  inc = zeros(nnos,1);  # Influence from incident plane waves

  for i=1:nnos # Loop over the source points
    xd=NOS[i,2]; # x coordinate of the source point
    yd=NOS[i,3]; # y coordinate of the source point
    for j=1:nelem # Loop over the elements
      noi::Int64=ELEM[j,2]; # First point of the element
      nof::Int64=ELEM[j,3]; # Second point of the element
      x1=NOS_GEO[noi,2]; # x coordinate of the first point of the element
      x2=NOS_GEO[nof,2]; # x coordinate of the second point of the element
      y1=NOS_GEO[noi,3]; # y coordinate of the first point of the element
      y2=NOS_GEO[nof,3];  # y coordinate of the second point of the element
      if i==j # The source point belongs to the element
        g,h = calcula_GeHs(x1,y1,x2,y2,k); 	# Singular integration
	qsitelles,Jtelles = telles(qsi,-1); # Evaluates the Telles' points and Jacobian
        g1, h1 = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k);	# Singular integration using the Telles transformation
	qsitelles,Jtelles = telles(qsi,1); # Evaluates the Telles' points and Jacobian
        g2, h1 = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k);	# Singular integration using the Telles transformation
	qsitelles,Jtelles = telles(qsi,0); # Evaluates the Telles' points and Jacobian
        g4, h4 = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k);	# Singular integration using the Telles transformation
	h1 = 0.5;
	g3 = g1 + g2;
        println("Diferença entre g e g3 = ", abs(g-g3))
        println("Diferença entre g e g4 = ", abs(g-g4))
      else # The source point doesnt belong to the element
        g,h = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsi,w,k);	# Non singular integration
      end
      G[i,j] = g
      H[i,j] = h
    end
  end
  # Rigid body consideration for the evaluation of the singular terms of matrix H
#for m = 1 : nnos
#    H[m,m] = 0; # let the diagonal be zero
#    for n = 1 : nnos	# Loop over the source points
#        if n != m	# If it's not a diagonal term
#            H[m,m] = H[m,m] - H[m,n]; # Subtracts the term from the diagonal
#        end;
#    end;
#end;

  return G,H
end

function cal_GeHpot(NOS,NOS_GEO,ELEM,k,fc,qsi,w)
# Evaluates the G and H matrices for the linear system H phi = G q, where phi is a vector containing the values of the velocity potential and q is a vector containing the values of the flux at the boundary of the problem.
 
  nelem::Int64=size(ELEM)[1]; # Number of elements
  nnos::Int64=nelem; # Number of nodes
  G=complex(zeros(nnos,nnos)); 	# Allocates matrix G
  H=complex(zeros(nnos,nnos));	# Allocates matrix H
  q=zeros(nnos,1);  # Influence from concentrated sources
  inc = zeros(nnos,1);  # Influence from incident plane waves
  qsitelles,Jtelles = telles(qsi,0); # Evaluates the Telles' points and Jacobian

  for i=1:nnos # Loop over the source points
    xd=NOS[i,2]; # x coordinate of the source point
    yd=NOS[i,3]; # y coordinate of the source point
    for j=1:nelem # Loop over the elements
      noi::Int64=ELEM[j,2]; # First point of the element
      nof::Int64=ELEM[j,3]; # Second point of the element
      x1=NOS_GEO[noi,2]; # x coordinate of the first point of the element
      x2=NOS_GEO[nof,2]; # x coordinate of the second point of the element
      y1=NOS_GEO[noi,3]; # y coordinate of the first point of the element
      y2=NOS_GEO[nof,3];  # y coordinate of the second point of the element
      if i==j # The source point belongs to the element
        #g,h = calcula_GeHs(x1,y1,x2,y2,k); 	# Singular integration
	g, h = calcula_GeHnspot(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k)
# Singular integration using the Telles transformation
	h = -0.5
        #htelles = htelles + 0.5	# Adding the jump term
        #println("Diferença entre g e gtelles = ", abs(g-gtelles))
        #println("Diferença entre h e htelles = ", abs(h-htelles))
      else # The source point doesnt belong to the element
        g,h = calcula_GeHnspot(x1,y1,x2,y2,xd,yd,qsi,w,k);	# Non singular integration
      end
      G[i,j] = g
      H[i,j] = h
    end
  end
  # Rigid body consideration for the evaluation of the singular terms of matrix H
#for m = 1 : nnos
#    H[m,m] = 0; # let the diagonal be zero
#    for n = 1 : nnos	# Loop over the source points
#        if n != m	# If it's not a diagonal term
#            H[m,m] = H[m,m] - H[m,n]; # Subtracts the term from the diagonal
#        end;
#    end;
#end;

  return G,H
end
