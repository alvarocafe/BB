function monta_phieq(CDC,x)
# Apply the boundary conditions on system A x = b, once the vector x has been determined, to separate the values of the velocity potential (phi) and flux (q), described as H phi = G q
  ncdc = length(CDC[:,1]);	# Number of boundary conditions
  nnos = length(x)	# Number of nodes
  phi = complex(zeros(nnos))	# Allocates the vector for the velocity potential
  q = complex(zeros(nnos))	# Allocates the vector for the flux
  for i=1:ncdc # Loop over the boundary conditions
    tipoCDC=CDC[i,2]; # Type of the boundary condition (Neumann or Dirichlet)
    valorCDC=CDC[i,3]; # Value of the boundary condition
    valorcalculado=x[i]; # Value which was previously unknown
    if tipoCDC == 1 # The flux is known
      phi[i] = valorcalculado; # The velocity potential was evaluated
      q[i] = valorCDC; # The flux is given by the boundary condition
    else # The velocity potential is known
      phi[i] = valorCDC; # The velocity potential is given by the boundary coindition
      q[i] = valorcalculado; # The flux was evaluated
    end
  end

  return phi,q
end
