function calc_fforma(qsi)
  # Calcula as funções de forma lineares contínuas N1 e N2
  N1=1/2*(1-qsi); # Função de forma N1 => linear contínua
  N2=1/2*(1+qsi); # Função de forma N2 => linear contínua
  return N1,N2
end

function calc_fforma_linear(qsi)
    #Calcula as funções de forma lineares descontínuas
    N1 = 2/4 - 3/4*qsi
    N2 = 2/4 + 3/4*qsi
    return [N1 N2]
end


function calc_fforma_quad(qsi)
  # Calcula as funções de forma quadráticas contínuas
  N1=(qsi/2)*(1-qsi); # Função de forma N1 => quadrática contínua
  N2=(1-qsi)*(1+qsi); # Função de forma N2 => quadrática contínua
  N3=(qsi/2)*(1+qsi); # Função de forma N3 => quadrática contínua
  N = [N1 N2 N3]
  return N
end

function calc_dfforma_quad(qsi)
  # Calcula as derivadas das funções de forma quadráticas contínua
  dN1=(2*qsi-1)/2; # Derivada da função de forma N1 => quadrática contínua
  dN2=-2*qsi; # Derivada da função de forma N2 => quadrática contínua
  dN3=(2*qsi +1)/2; # Derivada da função de forma N3 => quadrática contínua
  dN = [dN1 dN2 dN3]
  return dN
end

function calc_fforma_Bezier(qsi)
  # Calcula as funções de forma para o elemento de Bézier quadrático
  N1=(1-qsi)^2; # Função de forma N1 => quadrática contínua
  N2=2*qsi*(1-qsi); # Função de forma N2 => quadrática contínua
  N3=qsi^2; # Função de forma N3 => quadrática contínua
  N = [N1 N2 N3]
  return N
end

function calc_dfforma_Bezier(qsi)
  # Calcula as derivadas para as funções de forma para o elemento de Bézier quadrático
  dN1=2*(1-qsi); # Função de forma N1 => quadrática contínua
  dN2=2*(1-2*qsi); # Função de forma N2 => quadrática contínua
  dN3=2*qsi; # Função de forma N3 => quadrática contínua
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

function calc_fforma_triconst(qsi, eta)
  zeta = 1 - qsi - eta
  # Calcula as funções de forma triangulares lineares contínuas
  N1=qsi; # Função de forma N1
  N2=eta; # Função de forma N2
  N3=zeta; # Função de forma N3
  N = [N1 N2 N3]
  return N
end
function calc_dfforma_triconst(qsi, eta)
  # Calcula as derivadas para as funções de forma triangulares lineares contínuas
  dN1dqsi=1
  dN1deta=0
  dN2dqsi=0
  dN2deta=1
  dN3dqsi =-1
  dN3deta =-1
  dN = [dN1dqsi dN2dqsi dN3dqsi dN1deta dN2deta dN3deta]
  return dN
end
function cal_Jacobiano3D(x,y,z,dN)
    # Calcula o Jacobiano para elementos com 3 funções de forma 3D
    dxdqsi = dN[1]*x[1] + dN[2]*x[2] + dN[3]*x[3]
    dxdeta = dN[4]*x[1] + dN[5]*x[2] + dN[6]*x[3]
    dydqsi = dN[1]*y[1] + dN[2]*y[2] + dN[3]*y[3]
    dydeta = dN[4]*y[1] + dN[5]*y[2] + dN[6]*y[3]
    dzdqsi = dN[1]*z[1] + dN[2]*z[2] + dN[3]*z[3]
    dzdeta = dN[4]*z[1] + dN[5]*z[2] + dN[6]*z[3]
    J = sqrt((dydqsi*dzdeta - dzdqsi*dydeta)^2 + (dzdqsi*dxdeta - dxdqsi*dzdeta)^2 + (dxdqsi*dydeta - dydqsi*dxdeta)^2)
    return J
end

function calc_fforma_quadlin(qsi, eta)
  # Calcula as funções de forma quadrilaterais lineares contínuas
  N1=1/4*(1-qsi)*(1-eta) # Função de forma N1
  N2=1/4*(1+qsi)*(1-eta) # Função de forma N2
  N3=1/4*(1+qsi)*(1+eta) # Função de forma N3
  N4=1/4*(1-qsi)*(1+eta) # Função de forma N3
  N = [N1 N2 N3 N4]
  return N
end
function calc_dfforma_quadlin(qsi, eta)
  # Calcula as derivadas para as funções de forma triangulares lineares contínuas
  dN1dqsi=-1/4*(1-eta)
  dN2dqsi=1/4*(1-eta)
  dN3dqsi=1/4*(1+eta)
  dN4dqsi=-1/4*(1+eta)
  dN1deta=-1/4*(1-qsi)
  dN2deta=-1/4*(1+qsi)
  dN3deta=1/4*(1+qsi)
  dN4deta=1/4*(1-qsi)
  dN = [dN1dqsi dN2dqsi dN3dqsi dN4dqsi dN1deta dN2deta dN3deta dN4deta]
  return dN
end
function cal_Jacobiano3D_quadlin(x,y,z,dN)
    # Calcula o Jacobiano para elementos com 3 funções de forma 3D
    dxdqsi = dN[1]*x[1] + dN[2]*x[2] + dN[3]*x[3] + dN[4]*x[4]
    dxdeta = dN[5]*x[1] + dN[6]*x[2] + dN[7]*x[3] + dN[8]*x[4]
    dydqsi = dN[1]*y[1] + dN[2]*y[2] + dN[3]*y[3] + dN[4]*y[4]
    dydeta = dN[5]*y[1] + dN[6]*y[2] + dN[7]*y[3] + dN[8]*y[4]
    dzdqsi = dN[1]*z[1] + dN[2]*z[2] + dN[3]*z[3] + dN[4]*z[4]
    dzdeta = dN[5]*z[1] + dN[6]*z[2] + dN[7]*z[3] + dN[8]*z[4]
    J = sqrt((dydqsi*dzdeta - dzdqsi*dydeta)^2 + (dzdqsi*dxdeta - dxdqsi*dzdeta)^2 + (dxdqsi*dydeta - dydqsi*dxdeta)^2)
    return J
end

function calc_fforma_triquad(qsi, eta)
  zeta = 1 - qsi - eta
  # Calcula as funções de forma quadrilaterais lineares contínuas
  N1=qsi*(2*qsi -1) # Função de forma N1
  N2=eta*(2*eta -1) # Função de forma N2
  N3=zeta*(2*zeta -1) # Função de forma N3
  N4=4*qsi*eta # Função de forma N3
  N5=4*eta*zeta # Função de forma N3
  N6=4*qsi*zeta # Função de forma N3

  N = [N1 N2 N3 N4 N5 N6]
  return N
end
function calc_dfforma_triquad(qsi, eta)
  zeta = 1 - qsi - eta
  # Calcula as derivadas para as funções de forma triangulares lineares contínuas
  dN1dqsi=4*qsi-1
  dN2dqsi=0
  dN3dqsi=1-4*zeta
  dN4dqsi=4*eta
  dN5dqsi=-4*eta
  dN6dqsi=4*(zeta-qsi)

  dN1deta=0
  dN2deta=4*eta-1
  dN3deta=1-4*zeta
  dN4deta=4*qsi
  dN5deta=4*(zeta-eta)
  dN6deta=-4*qsi
  dN = [dN1dqsi dN2dqsi dN3dqsi dN4dqsi dN5dqsi dN6dqsi dN1deta dN2deta dN3deta dN4deta dN5deta dN6deta]
  return dN
end
function cal_Jacobiano3D_triquad(x,y,z,dN)
    # Calcula o Jacobiano para elementos com 3 funções de forma 3D
    dxdqsi = dN[1]*x[1] + dN[2]*x[2] + dN[3]*x[3] + dN[4]*x[4] + dN[5]*x[5] + dN[6]*x[6]
    dxdeta = dN[7]*x[1] + dN[8]*x[2] + dN[9]*x[3] + dN[10]*x[4]+ dN[11]*x[5] + dN[12]*x[6]
    dydqsi = dN[1]*y[1] + dN[2]*y[2] + dN[3]*y[3] + dN[4]*y[4] + dN[5]*y[5] + dN[6]*y[6]
    dydeta = dN[7]*y[1] + dN[8]*y[2] + dN[9]*y[3] + dN[10]*y[4]+ dN[11]*y[5] + dN[12]*y[6]
    dzdqsi = dN[1]*z[1] + dN[2]*z[2] + dN[3]*z[3] + dN[4]*z[4] + dN[5]*z[5] + dN[6]*z[6]
    dzdeta = dN[7]*z[1] + dN[8]*z[2] + dN[9]*z[3] + dN[10]*z[4]+ dN[11]*z[5] + dN[12]*z[6]
    J = sqrt((dydqsi*dzdeta - dzdqsi*dydeta)^2 + (dzdqsi*dxdeta - dxdqsi*dzdeta)^2 + (dxdqsi*dydeta - dydqsi*dxdeta)^2)
    return J
end

function calc_fforma_quadquad(qsi, eta)
  # Calcula as funções de forma quadrilaterais quadráticas contínuas
  N1=-1/4*(1-qsi)*(1-eta)*(qsi+1+eta) # Função de forma N1
  N2=1/4*(1+qsi)*(1-eta)*(qsi-1-eta) # Função de forma N2
  N3=1/4*(1+qsi)*(1+eta)*(eta-1+qsi) # Função de forma N3
  N4=1/4*(1-qsi)*(1+eta)*(-qsi-1+eta) # Função de forma N4

  N5=1/2*(1+qsi)*(1-qsi)*(1-eta) # Função de forma N5
  N6=1/2*(1+qsi)*(1+eta)*(1-eta) # Função de forma N6
  N7=1/2*(1+qsi)*(1+eta)*(-qsi+1) # Função de forma N7
  N8=1/2*(1-qsi)*(1+eta)*(1-eta) # Função de forma N8

  N = [N1 N2 N3 N4 N5 N6 N7 N8]
  return N
end
function calc_dfforma_quadquad(qsi, eta)
  # Calcula as derivadas para as funções de forma quadrilaterais quadráticas contínuas
  dN1dqsi=-1/4*(1-eta)
  dN2dqsi=1/4*(1-eta)
  dN3dqsi=1/4*(1+eta)
  dN4dqsi=-1/4*(1+eta)

  dN1deta=-1/4*(1-qsi)
  dN2deta=-1/4*(1+qsi)
  dN3deta=1/4*(1+qsi)
  dN4deta=1/4*(1-qsi)
  dN = [dN1dqsi dN2dqsi dN3dqsi dN4dqsi dN1deta dN2deta dN3deta dN4deta]
  return dN
end
function cal_Jacobiano3D_quadquad(x,y,z,dN)
    # Calcula o Jacobiano para elementos com 3 funções de forma 3D
    dxdqsi = dN[1]*x[1] + dN[2]*x[2] + dN[3]*x[3] + dN[4]*x[4]
    dxdeta = dN[5]*x[1] + dN[6]*x[2] + dN[7]*x[3] + dN[8]*x[4]
    dydqsi = dN[1]*y[1] + dN[2]*y[2] + dN[3]*y[3] + dN[4]*y[4]
    dydeta = dN[5]*y[1] + dN[6]*y[2] + dN[7]*y[3] + dN[8]*y[4]
    dzdqsi = dN[1]*z[1] + dN[2]*z[2] + dN[3]*z[3] + dN[4]*z[4]
    dzdeta = dN[5]*z[1] + dN[6]*z[2] + dN[7]*z[3] + dN[8]*z[4]
    J = sqrt((dydqsi*dzdeta - dzdqsi*dydeta)^2 + (dzdqsi*dxdeta - dxdqsi*dzdeta)^2 + (dxdqsi*dydeta - dydqsi*dxdeta)^2)
    return J
end
