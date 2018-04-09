#Input dos dados
p1 = [0 0]
p2 = [1 1]
p3 = [2 0]
NOS,ELEM,CDC = dad_quad(p1,p2,p3)

#Geração dos pontos e pesos de Gauss
npg=6*6;
qsi,w = Gauss_Legendre(0,1,npg)
n_pint = size(qsi,2)
N = calc_fforma_Bezier(0.5)
x_nocentral = N[1]*NOS[1,2]+N[2]*NOS[2,2]+N[3]*NOS[3,2]
y_nocentral = N[1]*NOS[1,3]+N[2]*NOS[2,3]+N[3]*NOS[3,3]
nodes = zeros(2,3)
nodes[1,:] = [NOS[1,2] x_nocentral NOS[3,2]]
nodes[2,:] = [NOS[1,3] y_nocentral NOS[3,3]]
x = [NOS[1,2] NOS[2,2] NOS[3,2]]
y = [NOS[1,3] NOS[2,3] NOS[3,3]]
#Integração dos elementos quadráticos
I1 = complex(0,0)
I2 = complex(0,0)
xx = zeros(size(qsi,2),1)
yy = zeros(size(qsi,2),1)
N = zeros(size(qsi,2),3)
dN = zeros(size(qsi,2),3)
dgamadqsi = zeros(size(qsi,2),1)
phiast =complex(zeros(size(qsi,2),1))
qast = complex(zeros(size(qsi,2),1))

for kk=1:n_pint # Laço sobre os pontos de integração
  N[kk,:] =calc_fforma_Bezier(qsi[kk]); # Calcula as funções de forma
  dN[kk,:] = calc_dfforma_Bezier(qsi[kk])
  dgamadqsi[kk]=cal_Jacobiano(x, y,dN[kk,:]); # Jacobiano
  xx[kk]=N[kk,1]*x[1]+N[kk,2]*x[2]+N[kk,3]*x[3]; # Calcula a coordenada x do ponto de integração
  yy[kk]=N[kk,1]*y[1]+N[kk,2]*y[2]+N[kk,3]*y[3]; # Calcula a coordenada y do ponto de integração
  dx=dN[kk,1]*x[1]+dN[kk,2]*x[2]+dN[kk,3]*x[3]; # Calcula a coordenada x do ponto de integração
  dy=dN[kk,1]*y[1]+dN[kk,2]*y[2]+dN[kk,3]*y[3]; # Calcula a coordenada y do ponto de integração
  sx=dx/dgamadqsi[kk]; # Componente x do vetor tangente
  sy=dy/dgamadqsi[kk]; # Componente y do vetor tangente
  nx=sy; # Componente x do vetor normal
  ny=-sx; # Componente y do vetor normal
  phiast[kk],qast[kk] =calc_solfund(xx[kk],yy[kk],nodes[1,1],nodes[2,1],nx,ny,340,340); # Obtemos as soluções fundamentais
  I1=I1+phiast[kk]*w[kk]*dgamadqsi[kk]; # Integral da matriz H
  I2=I2+qast[kk]*w[kk]*dgamadqsi[kk]; # Integral da matriz G
end
println("A integral I1= ",I1)
#Iniciamos a parte gráfica
close("all") #Vamos fechar as janelas que possam estar abertas

#Solução Fundamental
#figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
#plot(qsi',abs(phiast),label = "phiast",linestyle="-",linewidth=3,color ="green")
#plot(qsi',abs(qast),label = "qast",linestyle="-",linewidth=3,color="blue")
#legend(fontsize="12.0",loc="best")
#xlabel("xi")
#ylabel("Fundamental solutions")
#grid(1)
#title("Fundamental solution for the quadratic Bézier element")
#savefig("docs/discretização/figuras/fundamental_bezier.pdf",transparent = true)

#Funções de forma
#figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
#plot(qsi',N[:,1],label = "N1",linestyle="-",linewidth=3,color ="green")
#plot(qsi',N[:,2],label = "N2",linestyle="-",linewidth=3,color="blue")
#plot(qsi',N[:,3],label = "N3",linestyle="-",linewidth=3,color="red")
#legend(fontsize="12.0",loc="best")
#xlabel("xi")
#ylabel("Shape functions")
#grid(1)
#title("Shape functions for the quadratic Bézier element")
#savefig("docs/discretização/figuras/ff_formaquadratico_bezier.pdf",transparent = true)

#Jacobiano
#figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
#plot(qsi',dgamadqsi,label = "Jacobian",linestyle="-",linewidth=3,color ="red")
#legend(fontsize="12.0",loc="best")
#xlabel("xi")
#ylabel("Jacobian")
#grid(1)
#title("Jacobian for the quadratic Bézier element")
#savefig("docs/discretização/figuras/jacobian_bezier.pdf",transparent = true)


#Elemento
figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
plot(NOS[:,2],NOS[:,3],linestyle="-",linewidth=3,color="grey",marker="o",markersize=12,label="Geometrical points")
plot(xx,yy,marker=".",color="black",linestyle="--",markersize=12,label = "Gauss' points")
plot(nodes[1,:],nodes[2,:],marker="s",markersize=10,linestyle="none",color="blue",label = "Physical points (Nodes)")
grid(1)
PyPlot.xlabel("x",fontsize="12.0")
ylabel("y",fontsize="12.0")
title("Quadratic Bézier element",fontsize="16.0")
legend(fontsize="14.0",loc="lower center")
PyPlot.xlim(-0.2, 2.2)
PyPlot.ylim(-0.2, 1.2)
#savefig("docs/discretização/figuras/quadratico_Bezier.pdf",transparent = true)
