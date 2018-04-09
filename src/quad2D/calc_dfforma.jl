function calc_dfforma(qsi)
  # Calcula as derivadas das funções de forma quadráticas contínua
  dN1=(2*qsi-1)/2; # Derivada da função de forma N1 => quadrática contínua
  dN2=-2*qsi; # Derivada da função de forma N2 => quadrática contínua
  dN3=(2*qsi +1)/2; # Derivada da função de forma N3 => quadrática contínua
  dN = [dN1 dN2 dN3]
  return dN
end

