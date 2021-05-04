"""
Created on Thu Nov 19 11:18:59 2020

@author: fernando
"""
#caso (i)

using Plots, ColorSchemes
pyplot()

ne=200;
L=1;
k=kclosed(1,5);

phi_dom,PONTOS_dom=constclosed2D(ne,L,k)

n_pint = 100; # Number of domain points
pontosx = zeros(n_pint,2);
delta = 0.02; # distance from both ends 
passo = (L-2*delta)/(n_pint-1);
for i = 1:n_pint
    pontosx[i,:] = [i delta+i*passo];
end
pontosy = pontosx;

phi_real=real(phi_dom);
phi_contour = zeros(n_pint,n_pint);
for j = 1:n_pint
    index1=n_pint*(j-1)+1
    index2=n_pint*j
    phi_contour[:,j] = phi_real[index1:index2]
end

solar = ColorSchemes.solar.colors
contour(pontosx[:,2],pontosy[:,2],phi_contour,fill=true,colorbar=false,xaxis=("x"),yaxis=("y"),title=("Solução numérica de u caso (i)"),seriescolor=cgrad(ColorSchemes.viridis.colors),aspect_ratio=:equal)
