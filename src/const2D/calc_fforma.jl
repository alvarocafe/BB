function calc_fforma(qsi)
	# Evaluates the shape functions for continuous linear elements  
  N1=1/2*(1-qsi); # N1 - first shape function for the continuous linear element
  N2=1/2*(1+qsi); # N2 - second shape function for the continuous linear element
  return N1,N2
end
