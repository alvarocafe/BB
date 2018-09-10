function calc_fforma(qsi)
  # Evaluates the shape functions for the continuous quadratic element
  N1=(qsi/2)*(1-qsi); 
  N2=(1-qsi)*(1+qsi); 
  N3=(qsi/2)*(1+qsi); 
  N = [N1 N2 N3]
  return N
end

function calc_dfforma(qsi)
  # Calcula as derivadas das funções de forma quadráticas contínua
  dN1=(2*qsi-1)/2; # Derivada da função de forma N1 => quadrática contínua
  dN2=-2*qsi; # Derivada da função de forma N2 => quadrática contínua
  dN3=(2*qsi +1)/2; # Derivada da função de forma N3 => quadrática contínua
  dN = [dN1 dN2 dN3]
  return dN
end

function cal_Jacobiano(x,y,dN)
    # Calcula o Jacobiano para o elemento quadrático contínuo (Bézier ou Lagrangiano)
    dx = dN[1]*x[1] + dN[2]*x[2] + dN[3]*x[3]
    dy = dN[1]*y[1] + dN[2]*y[2] + dN[3]*y[3]
    J = sqrt(dx.^2 + dy.^2)
    return J
end

function telles(gamm,eet)

eest = eet^2 - 1;
term1 = eet*eest + abs(eest);
if term1 < 0
    term1 = (-term1)^(1/3);
    term1 = -term1;
else
    term1 = term1^(1/3);
end

term2 = eet*eest - abs(eest);
if term2 < 0
    term2 = (-term2)^(1/3);
    term2 = -term2;
else
    term2 = term2^(1/3);
end
GAMM = term1 + term2 + eet;


Q = 1 + 3*GAMM^2;
A = 1/Q;
B = -3*GAMM/Q;
C = 3*GAMM^2/Q;
D = -B;

eta = A*gamm.^3 + B*gamm.^2 + C*gamm + D;
Jt = 3*A*gamm.^2 + 2*B*gamm + C;
return  eta,Jt
end
