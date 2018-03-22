function aplica_CDC(G,H,CDC)
# Applies the boundary conditions by exchanging the columns of matrices H and G and constructs the linear system: Ax = b, where x contains all of the unknowns of the problem.
  ncdc = length(CDC[:,1]); # number of lines in the boundary conditions matrix
  A=H*1.0;	# Right side part of the linear system (unknowns)
  B=G*1.0;	# Left side part of the linear system (knowns)
  for i=1:ncdc # Loop on the number of boundary conditions
    tipoCDC = CDC[i,2]; # Type of the boundary condition
    if tipoCDC == 0 # The velocity potential is known
      colunaA=-A[:,i]; # 
      A[:,i]=-B[:,i]; # Matrix A receives the column from matrix G
      B[:,i]=colunaA; # Matrix b receives the column from matrix H
    end
  end
  valoresconhecidos=CDC[:,3]; # Boundary conditions' values
  b=B*valoresconhecidos; # array b
  return A, b
end


