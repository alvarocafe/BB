#Programa para resolver o problema do cilindro vibrante utilizando o BEM_base

# include("cal_GeH_linear.jl")
# include("docs/discretização/source/integra_linear.jl") # Exemplo da integração de um elemento linear

###Elemento constante################################
# Elemento constante,
# Input dos dados para reproduzir o problema utilizando elementos lineares descontínuos
#Geração dos pontos e pesos de Gauss
npg=16;
qsi,w = Gauss_Legendre(-1,1,npg)
# function compare_disc()
    include("dados/dad_2.jl")
    NOS_GEO,NOS,ELEM,CDC = format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg)
    nnos = size(NOS,1)
    b1 = 1:nnos
    # println("Construindo A e B")
    # A,B = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,CDC,1,1,qsi,w])
    # b = B*CDC[:,3]
    println("Construindo H e G")
    CW = 1
    FR = 1
    G,H = cal_HeG(NOS,NOS_GEO,ELEM,CW,FR,fc,qsi,w)
    A,b = aplica_CDC(G,H,CDC)
    x = A\b
    phi,qphi = monta_Teq(CDC,x)
#     println(phi)
#     println(NOS)
# return phi
# end
function cyl_const(i,FR)
  # Resolve o cilindro vibrante para elementos constantes
  # i = 20 # Numero de elementos por metade de círculo
  PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc,phi_analytical = dad_1(i,FR)
  # PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW, fc, finc,phi_analytical = dad_ext(i)
  NOS_GEO,NOS,ELEM,CDC = format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg)
  println("Número de elementos por  meio círculo = ",i)
  # Realizar a integração e resolver o sistema linear
  nnos = size(NOS,1)
  b1 = 1:nnos
  println("Construindo A e B")
  @time A,B = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,CDC,CW,FR,qsi,w])
  # println("Construindo H e G")
  # @time G,H = cal_HeG(NOS,NOS_GEO,ELEM,CW,FR,fc,qsi,w)
  # A,b = aplica_CDC(G,H,CDC)
  # println("Construindo b")
  b = B*CDC[:,3]
  #  println("Construindo x")
  x = A\b
  # println("Construindo phi, qphi")
  phi,qphi = monta_Teq(CDC,x)
  println("Construindo phi_pint")
  @time phi_pint = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,CW,FR,fc,finc,qsi,w)
  erro = abs((sum((phi_pint - phi_analytical).^2)))
  println("--------------------------------------------------------------------")
  println("--------------------------------------------------------------------")
  println("Coordenadas dos pontos de colocação:")
  println(NOS[:,2:3])
  #Parte gráfica
  # figure(figsize=(8,5))
  subplot("211")
  plot(NOS[:,2],NOS[:,3],linewidth=3,color="black",marker="o",markersize=3,label="Boundary elements")
  plot([NOS[1,2], NOS[nnos,2]],[NOS[1,3], NOS[nnos,3]],linewidth=3,color="black",marker="o",markersize=3)
  plot(PONTOS_int[:,2],PONTOS_int[:,3],linestyle="none",color="purple",marker="*",markersize=8,label="External points")
  ylabel("distance [m]",fontsize="18.0")
  grid("on")
  axis("equal")
  ax1 = gca()
  legend()
  subplot(212,sharex=ax1)
  plot(PONTOS_int[:,2],real.(phi_pint),linewidth=4,color="blue",marker="o",markersize=8,label="BEM")
  plot(PONTOS_int[:,2],real.(phi_analytical),linewidth=4,color="green",label="Analytical")
  legend()
  grid("on")
  xlabel("distance [m]",fontsize="18.0")
  ylabel("velocity potential []",fontsize="18.0")
  #savefig("docs/discretização/figuras/cilindros/cilindro_constante40.pdf",transparent = true)
  return erro,phi
end
# nerro = 20
# nfreq = 5
# erro = zeros(nerro,nfreq)
# x = zeros(nerro,nfreq)
# iter = linspace(100,1000,nerro)
# iterk = linspace(20,20000,nfreq)
# iq = 0
# for k in iterk
#     iq += 1
#     it = 0
#     for i in iter
#         it+=1
#         x[it,iq] = 2*i
#         erro[it,iq] = cyl_const(round(i),k)
#     end
# end
# #Parte gráfica
# close("all")
# #Erro
# figure(figsize=(6,5))
# loglog(x[:,1],erro[:,1],label="Frequency = 20 Hz")
# loglog(x[:,2],erro[:,2],label="Frequency = 5015 Hz")
# loglog(x[:,3],erro[:,3],label="Frequency = 10010 Hz")
# loglog(x[:,5],erro[:,4],label="Frequency = 15005 Hz")
# loglog(x[:,5],erro[:,5],label="Frequency = 20000 Hz")
# ylabel("Error")
# xlabel("Number of elements")
# #savefig("docs/discretização/figuras/erro_const2.pdf",transparent = true)
# #Elementos e soluções
# figure(figsize=(11,7))
# subplot("211")
# plot(NOS[:,2],NOS[:,3],linewidth=3,color="black",marker="o",markersize=3,label="Boundary elements")
# plot([NOS[1,2], NOS[nnos,2]],[NOS[1,3], NOS[nnos,3]],linewidth=3,color="black",marker="o",markersize=3)
# plot(PONTOS_int[:,2],PONTOS_int[:,3],linestyle="none",color="purple",marker="*",markersize=8,label="External points")
# ylabel("distance [m]",fontsize="18.0")
# grid("on")
# axis("equal")
# ax1 = gca()
# legend()
# subplot(212,sharex=ax1)
# plot(PONTOS_int[:,2],real.(phi_pint),linewidth=4,color="blue",marker="o",markersize=8,label="BEM")
# plot(PONTOS_int[:,2],real.(phi_analytical),linewidth=4,color="green",label="Analytical")
# legend()
# grid("on")
# xlabel("distance [m]",fontsize="18.0")
# ylabel("velocity potential []",fontsize="18.0")
# savefig("docs/discretização/figuras/cilindro_constante40.pdf",transparent = true)


# # #Elemento quadrático, resolvendo pelo cal_Aeb######################################################
# #Input dos dados
# NOS_GEO, NOS, ELEM, CDC, PONTOS_int, phi_analytical, CW, FR, fc, finc = dad_quad_cyl(4)
# ELEMconst = zeros(size(NOS,1),3);
# CDCconst = ones(size(NOS,1),3);
# for i = 1:size(NOS,1)
#   if i != size(NOS,1)
#     CDCconst[i,1] = i
#     ELEMconst[i,:] = [i i i+1]
#   else
#     CDCconst[i,1] = i
#     ELEMconst[i,:] = [i i 1]
#   end
#   CDCconst[:,3] = 1
# end
# #Geração dos pontos e pesos de Gauss
# npg=6;
# qsi,w = Gauss_Legendre(-1,1,npg)
# #Realizar a integração e resolver o sistema linear
# nnos = size(NOS,1)
# ne = size(ELEM,1)
# b1 = 1:nnos
# b2 = 1:ne
# A,B = cal_Aeb_quad(b1,b2, [NOS,NOS_GEO,ELEM,CDC,CW,FR,qsi,w])
# b = B*CDC[:,3]
# x = A\b
# phi,qphi = monta_Teq(CDC,x)
# phi_pint = calc_phi_pint_quad(PONTOS_int,NOS_GEO,ELEM,phi,qphi,CW,FR,fc,finc,qsi,w)
#
# # Para elementos constantes com os mesmos nos
# Aconst,Bconst = cal_Aeb(b1,b1, [NOS,NOS,ELEMconst,CDCconst,CW,FR,qsi,w])
# bconst = Bconst*CDCconst[:,3]
# xconst = Aconst\bconst
# phiconst,qphiconst = monta_Teq(CDCconst,xconst)
# phi_pintconst = calc_phi_pint(PONTOS_int,NOS_GEO,ELEMconst,phiconst,qphiconst,CW,FR,fc,finc,qsi,w)
#
# #Parte gráfica
# # close("all")
# figure(figsize=(11,7))
# title("Vibrating cylinder for quadratic elements")
# subplot("211")
# plot(NOS[:,2],NOS[:,3],linewidth=3,color="black",marker="o",markersize=3)
# ylabel("distance [m]",fontsize="18.0")
# grid("on")
# axis("equal")
# ax1 = gca()
# subplot(212,sharex=ax1)
# plot(PONTOS_int[:,2],abs.(phi_pint),linewidth=4,color="blue",marker="o",markersize=8,label="BEM")
# plot(PONTOS_int[:,2],abs.(phi_analytical),linewidth=4,color="green",label="Analytical")
# legend()
# grid("on")
# xlabel("distance [m]",fontsize="18.0")
# ylabel("velocity potential []",fontsize="18.0")
# # savefig("docs/cilindro vibrante/figuras/cilindro_quadratico.pdf",transparent = true)


#########################################################
# #Elemento de Bezier quadrático, resolvendo pelo cal_Aeb
# #Input dos dados
# p1 = [0 0]
# p2 = [1 1]
# p3 = [2 0]
# NOS,ELEM,CDC = dad_quad(p1,p2,p3)
# NOS_GEO = NOS
# CW = 1
# FR = 1
# #Geração dos pontos e pesos de Gauss
# npg=6*6;
# qsi,w = Gauss_Legendre(-1,1,npg)
# # Realizar a integração e resolver o sistema linear
# nnos = size(NOS,1)
# b1 = 1:nnos
# A,B = cal_Aeb_Bezier(b1,1, [NOS,NOS_GEO,ELEM,CDC,CW,FR,qsi,w])
# b = B*CDC[:,3]
# x = A\b
# phi,qphi = monta_Teq(CDC,x)

########################################################
# #Elemento NURBS, resolvendo pelo cal_HeG
#Input dos dados
function cyl_nurbs(h,FR)
PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc,phi_analytical = dad_1(1,FR)
crv = format_dad_iso(PONTOS,SEGMENTOS,MALHA)
dcrv=map(x->nrbderiv(x),crv)
n = length(crv);	# N�mero total de elementos
# h=10;#refinamento h
if h != 0
for i=1:n
  novosnos=linspace(0,1,h+2)
  degree=crv[i].order-1
  coefs,knots = bspkntins(degree,crv[i].coefs,crv[i].knots,novosnos[2:end-1])
  crv[i] = nrbmak(coefs,knots)
end
end
z=0;
for k=1:n
  for i=1:crv[k].number
    z=z+1
  end
end
numcurva=zeros(Integer,z)
collocPts=zeros(z)
CDC=zeros(z,3)
collocCoord=zeros(z,3)
z=0;
nnos=zeros(Integer,n)
for k=1:n
  p=crv[k].order-1;
  nnos[k]=crv[k].number;
  valorCDC=CCSeg[k,3];
  tipoCDC=CCSeg[k,2];
  for i=1:crv[k].number
    z=z+1;
    numcurva[z]=k;
    collocPts[z]=sum(crv[k].knots[(i+1):(i+p)])/p;  # Determinar a posição no vetor knots onde será feita a colocação
    if(i==2)
      collocPts[z-1]=(collocPts[z]+collocPts[z-1])/2;
    end
    if(i==nnos[k])
      collocPts[z]=(collocPts[z]+collocPts[z-1])/2;
    end
    CDC[z,:] = [z,tipoCDC,valorCDC];
  end
end
nnos2=cumsum([0 nnos'],2);

E=zeros(length(collocPts),length(collocPts));
for i=1:length(collocPts)
  collocCoord[i,:]=nrbeval(crv[numcurva[i]], collocPts[i]);
  B, id = nrbbasisfun(crv[numcurva[i]],collocPts[i])
  E[i,id+nnos2[numcurva[i]]]=B
end
println("Coordenadas dos pontos de colocação: ")
println(collocCoord)
H, G = calcula_iso(collocCoord,nnos2,crv,dcrv,FR/CW)
H=H+E/2;
#println(sum(H[9,:])) # Mostra a soma da linha para o elemento de comparação
A,b= aplica_CDCiso(G,H,CDC,E);
x=A\b; # Calcula o vetor x
phic,qc=monta_Teqiso(CDC,x); # Separa temperatura e fluxo
phi=E*phic;
q=E*qc;
# # Calculando o potencial de velocidade nos pontos internos
npg=6;
qsi,w = Gauss_Legendre(-1,1,npg)
phi_pint = calc_phi_pint_nurbs(collocCoord,nnos2,PONTOS_int,crv,dcrv,FR/CW,qsi,w,fc,finc,phi,q)
#Parte gráfica
# close("all")
# figure(figsize=(8,5))
title("Vibrating cylinder for NURBS elements")
subplot("211")
mostra_geo(crv)
plot(collocCoord[:,1],collocCoord[:,2],marker="s",markersize=10,linestyle="none",color="blue",label = "Physical points (Nodes)")
legend(fontsize="14.0",loc="best")
ylabel("distance [m]",fontsize="18.0")
grid("on")
axis("equal")
ax1 = gca()
PyPlot.xlim(-0.6, 0.6)
PyPlot.ylim(-0.6, 0.6)
subplot(212,sharex=ax1)
# subplot(212,sharex=ax1)
plot(PONTOS_int[:,2],abs.(phi_pint),linewidth=4,color="blue",marker="o",markersize=8,label="BEM")
plot(PONTOS_int[:,2],abs.(phi_analytical),linewidth=4,color="green",label="Analytical")
legend()
grid("on")
xlabel("distance [m]",fontsize="18.0")
ylabel("velocity potential []",fontsize="18.0")
# savefig("docs/discretização/figuras/cilindro_nurbs.pdf",transparent = true)
return phi
end
