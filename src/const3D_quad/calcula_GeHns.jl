function calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xd,yd,zd,n,qsi,w,k)
#integra��o n�o singular
n_pint=length(qsi); # N�mero de pontos de integra��o.
g=complex(0,0); # Inicializa o somatorio de g
h=complex(0,0); # Inicializa o somatorio de h
N = zeros(4)

for l=1:n_pint # La�o sobre os pontos de integra��o
    for m=1:n_pint # La�o sobre os pontos de integra��o
#	N=calc_fforma(qsi[l],qsi[m]); #  fun��es de forma
	N = (1./4.)*[(1. - qsi[l]).*(1. - qsi[m]);
               (1. + qsi[l]).*(1. - qsi[m]);
               (1. + qsi[l]).*(1. + qsi[m]);
               (1. - qsi[l]).*(1. + qsi[m])];
        x=N[1]*x1+N[2]*x2+N[3]*x3+N[4]*x4; # coordenada x do ponto de integra��o
        y=N[1]*y1+N[2]*y2+N[3]*y3+N[4]*y4; # coordenada y do ponto de integra��o
        z=N[1]*z1+N[2]*z2+N[3]*z3+N[4]*z4; # coordenada z do ponto de integra��o
        #J = calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,qsi[l],qsi[m]);# jacobiano
	#dNdqsi, dNdeta = calc_dfforma(qsi[l],qsi[m]); # Calcula a derivada das fun��es de forma
	dNdqsi = (1./4.)*[-(1. - qsi[m]);
		          (1. - qsi[m]);
		          (1. + qsi[m]);
		         -(1. + qsi[m])];

	dNdeta = (1./4.)*[-(1. - qsi[l]);
		         -(1. + qsi[l]);
		          (1. + qsi[l]);
		          (1. - qsi[l])];

	dxdqsi = x1*dNdqsi[1]+x2*dNdqsi[2]+x3*dNdqsi[3]+x4*dNdqsi[4];
	dydqsi = y1*dNdqsi[1]+y2*dNdqsi[2]+y3*dNdqsi[3]+y4*dNdqsi[4];
	dzdqsi = z1*dNdqsi[1]+z2*dNdqsi[2]+z3*dNdqsi[3]+z4*dNdqsi[4];

	dxdeta = x1*dNdeta[1]+x2*dNdeta[2]+x3*dNdeta[3]+x4*dNdeta[4];
	dydeta = y1*dNdeta[1]+y2*dNdeta[2]+y3*dNdeta[3]+y4*dNdeta[4];
	dzdeta = z1*dNdeta[1]+z2*dNdeta[2]+z3*dNdeta[3]+z4*dNdeta[4];

	g1 = dydqsi*dzdeta - dzdqsi*dydeta;
	g2 = dzdqsi*dxdeta - dxdqsi*dzdeta;
	g3 = dxdqsi*dydeta - dydqsi*dxdeta;
	J = sqrt(g1^2.0 + g2^2.0 + g3^2.0);

        #Tast,qast=calc_solfund(x,y,z,xd,yd,zd,n,k); # Solu��es fundamentais
	rx=x-xd;
	ry=y-yd;
	rz=z-zd;

	r =sqrt(rx^2+ry^2+rz^2);
	ZW=complex(0.,-k*r);
	Tast=exp(ZW)/complex(4*pi*r,0.);
	drdn=(rx*n[1] + ry*n[2] + rz*n[3])/r;
	qast=(ZW-complex(1.,0.))*Tast*complex(drdn/r,0.);

#	Tast = 1 #To test the numerical integration, the value of the fundamental solutions is set to 1. In this special case, the value of the sum of the rows of matrices H or G should equal the surface area of the model.
#	qast = 1        
	g=g+Tast*complex(w[l]*w[m]*J,0); # Integral da matriz G
        h=h+qast*complex(w[l]*w[m]*J,0); # Integral da matriz H
    end
end
return g,h
end

