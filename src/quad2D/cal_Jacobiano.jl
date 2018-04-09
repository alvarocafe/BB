function cal_Jacobiano(x,y,dN)
    # Calcula o Jacobiano para o elemento quadrático contínuo (Bézier ou Lagrangiano)
    dx = dN[1]*x[1] + dN[2]*x[2] + dN[3]*x[3]
    dy = dN[1]*y[1] + dN[2]*y[2] + dN[3]*y[3]
    J = sqrt(dx.^2 + dy.^2)
    return J
end

