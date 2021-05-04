"""
Created on Wed Nov 18 11:20:11 2020

@author: fernando
"""

#caso (i)

using Plots
pyplot()

ne=50
L=1
k=kclosed(1,5)

phi_dom,PONTOS_dom=constclosed2D(ne,L,k)

plot(PONTOS_dom[:,2],PONTOS_dom[:,3],real(phi_dom),st=:surface,c=:viridis,xaxis=("x"),yaxis=("y"),title=("Solução numérica de u caso (i)"))

#plot(PONTOS_dom[:,2],PONTOS_dom[:,3],phi_closed(kclosed(1,5),PONTOS_dom[:,2]),st=:surface,c=:viridis,xaxis=("x"),yaxis=("y"),title=("Solução analítica de \$\\phi\$ caso (i)"))

