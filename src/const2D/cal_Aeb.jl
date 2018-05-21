function cal_Aeb(b1,b2,arg)
# Builds the matrices for the linear system A x = b 
  NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k=arg
  nelem::Int64=size(ELEM)[1]; # Number of elements
  qsitelles,Jtelles = telles(qsi,0); # Evaluating the Telles' points and its Jacobian
  B=complex(zeros(length(b1),length(b2)));	# Allocates matrix B
  A=complex(zeros(length(b1),length(b2)));	# Allocates matrix A
  q=zeros(length(b1),1);	# Allocates array q
  ci=0
  for i in b1 # Loop over the source points
    ci+=1
    xd=NOS[i,2]; # x coordinate of the source point
    yd=NOS[i,3]; # y coordinate of the source point
    cj=0
    for j in b2 # Loop over the elements
      cj+=1
      noi::Int64=ELEM[j,2]; # Initial point of the element
      nof::Int64=ELEM[j,3]; # Final point of the element
      x1=NOS_GEO[noi,2]; # x coordinate of the first point of the element
      x2=NOS_GEO[nof,2]; # x coordinate of the second point of the element
      y1=NOS_GEO[noi,3]; # y coordinate of the first point of the element
      y2=NOS_GEO[nof,3];  # y coordinate of the second point of the element
      if i==j # The source point belongs to the element
        #g,h = calcula_GeHs(x1,y1,x2,y2,1,k);	# Singular integration
	g,h = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k);	# Singular integration using the Telles transformation
	h = 0.5;
        # println("Diferença entre g e gtelles = ", abs(g-gtelles))
        # println("Diferença entre h e htelles = ", abs(h-htelles))
      else # O ponto fonte n�o pertence ao elemento
        g,h = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsi,w,k); # Non singular integration
      end
	# Applying the boundary conditions
      if CDC[j,2]==0	# The velocity potential is known
        B[ci,cj] = -h	# Matrix B receives the value from matrix h
        A[ci,cj] = -g	# Matrix A receives the value from matrix g
      else
        B[ci,cj] = g	# Matrix B receives the value from matrix g
        A[ci,cj] = h	# Matrix A receives the value from matrix h
      end
    end
  end
	b = B*(CDC[b2,3])  # Builds the b array for the linear system
return A,b
end

function cal_Aeb_POT(b1,b2,arg)
# Builds the matrices for the linear system A x = b 
  NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k=arg
  nelem::Int64=size(ELEM)[1]; # Number of elements
  qsitelles,Jtelles = telles(qsi,0); # Evaluating the Telles' points and its Jacobian
  B=complex(zeros(length(b1),length(b2)));	# Allocates matrix B
  A=complex(zeros(length(b1),length(b2)));	# Allocates matrix A
  q=zeros(length(b1),1);	# Allocates array q
  ci=0
  for i in b1 # Loop over the source points
    ci+=1
    xd=NOS[i,2]; # x coordinate of the source point
    yd=NOS[i,3]; # y coordinate of the source point
    cj=0
    for j in b2 # Loop over the elements
      cj+=1
      noi::Int64=ELEM[j,2]; # Initial point of the element
      nof::Int64=ELEM[j,3]; # Final point of the element
      x1=NOS_GEO[noi,2]; # x coordinate of the first point of the element
      x2=NOS_GEO[nof,2]; # x coordinate of the second point of the element
      y1=NOS_GEO[noi,3]; # y coordinate of the first point of the element
      y2=NOS_GEO[nof,3];  # y coordinate of the second point of the element
      if i==j # The source point belongs to the element
        #g,h = calcula_GeHs(x1,y1,x2,y2,1,k);	# Singular integration
	g,h = calcula_GeHnspot(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k);	# Singular integration using the Telles transformation
	h = -0.5;
        # println("Diferença entre g e gtelles = ", abs(g-gtelles))
        # println("Diferença entre h e htelles = ", abs(h-htelles))
      else
        g,h = calcula_GeHnspot(x1,y1,x2,y2,xd,yd,qsi,w,k); # Non singular integration
      end
	# Applying the boundary conditions
      if CDC[j,2]==0	# The temperature is known
        B[ci,cj] = -h	# Matrix B receives the value from matrix h
        A[ci,cj] = -g	# Matrix A receives the value from matrix g
      else
        B[ci,cj] = g	# Matrix B receives the value from matrix g
        A[ci,cj] = h	# Matrix A receives the value from matrix h
      end
    end
  end
	b = B*(CDC[b2,3])  # Builds the b array for the linear system
return A,b
end
