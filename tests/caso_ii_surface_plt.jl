"""
Created on Thu Nov 19 10:46:49 2020

@author: fernando
"""
#caso (ii)

using Plots
pyplot()

ne=50
L=1
n=5
k=(2*n+1)*π./(2*L)

phi_dom,PONTOS_dom=cup2D(ne,L,k)

plot(PONTOS_dom[:,2],PONTOS_dom[:,3],real(phi_dom),st=:surface,c=:viridis,xaxis=("x"),yaxis=("y"),title=("Solução numérica de u caso (ii)"))

#plot(PONTOS_dom[:,2],PONTOS_dom[:,3],phi_cup(kcup(1,5),PONTOS_dom[:,2]),st=:surface,c=:viridis,xaxis=("x"),yaxis=("y"),title=("Solução analítica de u caso (ii)"))
