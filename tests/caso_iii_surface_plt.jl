#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 11:02:40 2020

@author: fernando
"""

#caso (iii)

using Plots
pyplot()

ne=50
L=1
k=kclosed(1,2)

phi_dom,PONTOS_dom=open2D(ne,L,k)

max_phi_dom=max(real(phi_dom)...)
normalized_phi_dom=real(phi_dom)./(max_phi_dom)

plot(PONTOS_dom[:,2],PONTOS_dom[:,3],normalized_phi_dom,st=:surface,c=:viridis,xaxis=("x"),yaxis=("y"),title=("Solução numérica de u caso (iii)"))

#plot(PONTOS_dom[:,2],PONTOS_dom[:,3],phi_closed(kclosed(1,7),PONTOS_dom[:,2]),st=:surface,c=:viridis,xaxis=("x"),yaxis=("y"),title=("Solução analítica de \$\\phi\$ caso (i)"))
