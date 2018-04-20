function calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xd,yd,zd,n,qsi,w,k)
#integra��o n�o singular
n_pint=length(qsi); # N�mero de pontos de integra��o.
g=complex(0,0); # Inicializa o somatorio de g
h=complex(0,0); # Inicializa o somatorio de h

tff = 0.0;
tj = 0.0;
ts = 0.0;
for l=1:n_pint # La�o sobre os pontos de integra��o
    for m=1:n_pint # La�o sobre os pontos de integra��o
	tic()
	N1,N2,N3,N4=calc_fforma(qsi[l],qsi[m]); #  fun��es de forma
	tff+=toq();
        x=N1*x1+N2*x2+N3*x3+N4*x4; # coordenada x do ponto de integra��o
        y=N1*y1+N2*y2+N3*y3+N4*y4; # coordenada y do ponto de integra��o
        z=N1*z1+N2*z2+N3*z3+N4*z4; # coordenada z do ponto de integra��o
	tic()
        J = calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,qsi[l],qsi[m]);# jacobiano
	tj+=toq();
	tic()
        Tast,qast=calc_solfund(x,y,z,xd,yd,zd,n,k); # Solu��es
        #  fundamentais
	ts+=toq()
#	Tast = 1 #To test the numerical integration, the value of the fundamental solutions is set to 1. In this special case, the value of the sum of the rows of matrices H or G should equal the surface area of the model.
#	qast = 1        
	g=g+Tast*complex(w[l]*w[m]*J,0); # Integral da matriz G
        h=h+qast*complex(w[l]*w[m]*J,0); # Integral da matriz H
    end
end
return g,h,tff,tj,ts
end

