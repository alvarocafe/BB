function calc_fforma(qsi)
  # Evaluates the shape functions for the continuous quadratic element
  N1=(qsi/2)*(1-qsi); 
  N2=(1-qsi)*(1+qsi); 
  N3=(qsi/2)*(1+qsi); 
  N = [N1 N2 N3]
  return N
end

