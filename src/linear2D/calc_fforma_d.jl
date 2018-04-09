function calc_fforma_d(qsi)
# Evaluates the discontinuous linear shape functions for points located at a quarter from the endpoints of the element.
N1=1/2 -1*qsi;
N2=1/2 + 1*qsi;
return N1,N2
end
