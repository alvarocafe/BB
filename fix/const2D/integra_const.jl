#Input dos dados
p1 = [0 0]
p2 = [1 1]
NOS_GEO,NOS,ELEM,CDC = dad_const(p1,p2)
CW = 1
FR = 1

#Geração dos pontos e pesos de Gauss
npg=6*6;
qsi,w = Gauss_Legendre(-1,1,npg)
qsitelles,Jtelles = telles(qsi,1); # Calculando os pontos de Telles e o Jacobiano

#Integração dos elementos constantes
nelem = size(ELEM,1)
nnos = size(NOS,1)
I1 = complex(0,0)
I2 = complex(0,0)
I1telles = complex(0,0)
I2telles = complex(0,0)
x = zeros(npg)
y = zeros(npg)
N1 = zeros(npg)
N2 = zeros(npg)
dgamadqsi = zeros(npg)
phiast =complex(zeros(npg))
qast = complex(zeros(npg))
phiasttelles =complex(zeros(npg))
qasttelles = complex(zeros(npg))

for j = 1:nelem

  no_inicial = ELEM[j,2]
  no_final = ELEM[j,3]

  for i = 1:npg
    L=sqrt((NOS_GEO[no_final,2]-NOS_GEO[no_inicial,2])^2+(NOS_GEO[no_final,3]-NOS_GEO[no_inicial,3])^2); # Comprimento do elemento
    dgamadqsi[i]=L/2; # Jacobiano
    sx=(NOS_GEO[no_final,2]-NOS_GEO[no_inicial,2])/L; # Componente x do vetor tangente
    sy=(NOS_GEO[no_final,3]-NOS_GEO[no_inicial,3])/L; # Componente y do vetor tangente
    nx=sy; # Componente x do vetor normal
    ny=-sx; # Componente y do vetor normal
    #Comparar Gauss com Telles
    N1[i],N2[i] = calc_fforma(qsi[i])
    x[i] = NOS_GEO[no_inicial,2]*N1[i] + NOS_GEO[no_final,2]*N2[i]
    y[i] = NOS_GEO[no_inicial,3]*N1[i] + NOS_GEO[no_final,3]*N2[i]
    phiast[i],qast[i] = calc_solfund(x[i],y[i],NOS_GEO[2,2],NOS_GEO[2,3],nx,ny,1,1)

    N1[i],N2[i] = calc_fforma(qsitelles[i])
    x[i] = NOS_GEO[no_inicial,2]*N1[i] + NOS_GEO[no_final,2]*N2[i]
    y[i] = NOS_GEO[no_inicial,3]*N1[i] + NOS_GEO[no_final,3]*N2[i]
    phiasttelles[i],qasttelles[i] = calc_solfund(x[i],y[i],NOS_GEO[2,2],NOS_GEO[2,3],nx,ny,1,1)
    # phiast[i],qast[i] = [1 1]
    I1telles = I1telles + w[i]*phiasttelles[i]*dgamadqsi[i]*Jtelles[i]
    I2telles = I2telles + w[i]*qasttelles[i]*dgamadqsi[i]*Jtelles[i]
    I1 = I1 + w[i]*phiast[i]*dgamadqsi[i]
    I2 = I2 + w[i]*qast[i]*dgamadqsi[i]
  end
end
println("O valor das integrais I1 e I2 são: I1= ",I1," e I2= ",I2)
println("O valor das integrais I1telles e I2telles são: I1telles= ",I1telles," e I2telles= ",I2telles)

fc = zeros(1,1)

G,H = cal_HeG(NOS,NOS_GEO,ELEM,CW,FR,fc,qsi,w)
A,b = aplica_CDC(G,H,CDC)
b1 = 1:nnos; b2 = 1:nelem;
A1,b1 = cal_Aeb(b1,b2,[NOS,NOS_GEO,ELEM,CDC,CW,FR,qsi,w])
println("Diferença entre A e A1 = ", abs(A-A1))
println("Diferença entre b e b1 = ", abs(b-b1))

#Iniciamos a parte gráfica
close("all") #Vamos fechar as janelas que possam estar abertas

#Solução Fundamental
figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
plot(qsi,abs(phiast),"ro",linestyle="-",linewidth=2,label = L" $\phi^{*}$ Gauss")
plot(qsi,abs(qast),"ro",linestyle="-",linewidth=2,label = L"$q^{*}$ Gauss")
plot(qsitelles,abs(phiasttelles),"-b>",markersize=6,linewidth=1,label = L" $\phi^{*}$ Telles")
plot(qsitelles,abs(qasttelles),"-b>",markersize=6,linewidth=1,label = L" $q^{*}$ Telles")
legend(fontsize="12.0",loc="best")
xlabel(L"\xi")
ylabel("Fundamental solutions")
grid(1)
title("Absolute value of the fundamental solution for the linear element")
# PyPlot.xlim(-0.2, 1.2)
PyPlot.ylim(-0.1, 1.65)
savefig("docs/discretização/figuras/fundamental_linear_telles2.pdf",transparent = true)


#Funções de forma
#figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
#plot(qsi',N1,label = "N1",linestyle="-",linewidth=3,color ="green")
#plot(qsi',N2,label = "N2",linestyle="-",linewidth=3,color="blue")
#legend(fontsize="12.0",loc="best")
#xlabel("xi")
#ylabel("Shape functions")
#grid(1)
#title("Shape functions for the linear element")
#savefig("docs/discretização/figuras/ff_forma_const.pdf",transparent = true)

#Jacobiano
#figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
#plot(qsi',dgamadqsi,label = "Jacobian",linestyle="-",linewidth=3,color ="red")
#legend(fontsize="12.0",loc="best")
#xlabel("xi")
#ylabel("Jacobian")
#grid(1)
#title("Jacobian for the constant element")
#savefig("docs/discretização/figuras/jacobian_const.pdf",transparent = true)


#Elemento
figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
plot(NOS_GEO[:,2],NOS_GEO[:,3],linestyle="-",linewidth=3,color="grey",marker="o",markersize=12,label="Geometrical points")
plot(x,y,marker=".",color="black",linestyle="--",markersize=12,label = "Telles' points")
plot(NOS_GEO[:,2],NOS_GEO[:,3],marker="s",markersize=10,linestyle="none",color="blue",label = "Physical points (Nodes)")
grid(1)
PyPlot.xlabel("x",fontsize="12.0")
ylabel("y",fontsize="12.0")
title("Linear element",fontsize="16.0")
legend(fontsize="14.0",loc="best")
PyPlot.xlim(-0.2, 1.2)
PyPlot.ylim(-0.2, 1.2)
savefig("docs/discretização/figuras/linearTelles2.pdf",transparent = true)
