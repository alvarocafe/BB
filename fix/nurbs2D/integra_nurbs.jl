#Input dos dados
p1 = [0 0]
p2 = [1 0]
println("Integrando uma curva NURBS.")
collocCoord,nnos,crv,dcrv,CDC,E = dad_iso(10)
H, G = calcula_iso(collocCoord,nnos,crv,dcrv,1)
H=H+E/2;
#println(sum(H[9,:])) # Mostra a soma da linha para o elemento de comparação
A,b= aplica_CDCiso(G,H,CDC,E);
x=A\b; # Calcula o vetor x
Tc,qc=monta_Teqiso(CDC,x); # Separa temperatura e fluxo
T=E*Tc;
q=E*qc;
close("all")
mostra_geo(crv)
plot(collocCoord[:,1],collocCoord[:,2],marker="s",markersize=10,linestyle="none",color="blue",label = "Physical points (Nodes)")
axis("equal")
grid(1)
PyPlot.xlabel("x",fontsize="12.0")
PyPlot.ylabel("y",fontsize="12.0")
title("NURBS element",fontsize="16.0")
legend(fontsize="14.0",loc="best")
PyPlot.xlim(-0.1, 1.1)
PyPlot.ylim(-0.6, 0.1)
savefig("docs/discretização/figuras/nurbs1.pdf",transparent = true)
