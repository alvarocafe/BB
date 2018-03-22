function trocacol(G,H,CDC)
# Applies the boundary conditions by exchanging the columns of matrices H and G
  ncdc = length(CDC); # Number of boundary conditions
  A=H*1.0;	# Make a copy of H
  B=G*1.0;	# Make a copy of G
  for i=1:ncdc # Loop over the boundary conditions
    tipoCDC = CDC[i]; # Type of boundary condition
    if tipoCDC == 0 # The velocity potential is known
      colunaA=-A[:,i]; # 
      A[:,i]=-B[:,i]; # Matrix A receives the column from matrix G
      B[:,i]=colunaA; # Matrix B receives the column from matrix H
    end
  end
#  b=B*valoresconhecidos; # Vector b
  return A, B
end
