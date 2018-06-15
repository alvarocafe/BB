function mostra_problema(ELEM,NOS_GEO,NOS,tipoCDC,valorCDC,normal)
nelem=size(ELEM,1);
figure();
ind=collect(1:nelem);
indKnownq=ind[tipoCDC];
indKnownT=ind[!tipoCDC];
maxT=maximum(abs(T[indKnownT]))
maxq=maximum(abs(q[indKnownq]))
xmax=maximum(NOS_GEO[:,1])
ymax=maximum(NOS_GEO[:,2])
xmin=minimum(NOS_GEO[:,1])
ymin=minimum(NOS_GEO[:,2])
deltax=xmax-xmin
deltay=ymax-ymin
dmax=sqrt(deltax^2+deltay^2)
ax = plt.gca() # get current axes
plt.grid("on")
plt.quiver(NOS[indKnownT,1],NOS[indKnownT,2],normal[1,indKnownT].*valorCDC[indKnownT],
      normal[2,indKnownT].*valorCDC[indKnownT],color="red",width=0.002,scale=50, headaxislength=0)
plt.quiver(NOS[indKnownq,1],NOS[indKnownq,2],normal[1,indKnownq].*valorCDC[indKnownq],
      normal[2,indKnownq].*valorCDC[indKnownq],color="blue",width=0.002,scale=50)
plt.plot(NOS[indKnownT,1],NOS[indKnownT,2],"ro",markersize=4);	# Plot the node of the elements
plt.plot(NOS[indKnownq,1],NOS[indKnownq,2],"bo",markersize=4);	# Plot the node of the elements
ELEM2=[ELEM ELEM[:,2]];
plt.triplot(NOS_GEO[:,1], NOS_GEO[:,2], ELEM2-1, color=(0.0,0.,0.),linewidth=0.4) 	# Plot the countour of the problem
plt.axis("equal")
ax[:set_xlim]((xmin-.15*deltax,xmax+.15*deltax));
ax[:set_ylim]((ymin-.15*deltay,ymax+.15*deltay));
end

function mostra_heatmap(NOS,PONTOS_INT,T,Ti,NOS_GEO,ELEM,dTdx,dTdy)
nelem=size(ELEM,1);
figure();
xmax=maximum(NOS_GEO[:,1])
ymax=maximum(NOS_GEO[:,2])
xmin=minimum(NOS_GEO[:,1])
ymin=minimum(NOS_GEO[:,2])
deltax=xmax-xmin
deltay=ymax-ymin
dmax=sqrt(deltax^2+deltay^2)

XY=[NOS;PONTOS_INT];
triang = tri.Triangulation(XY[:,1], XY[:,2])
cor=[T; Ti]
t=triang[:triangles]+1
centroid=(XY[t[:,1],:]+XY[t[:,2],:]+XY[t[:,3],:])/3.0;
i = inpoly(centroid,NOS,ELEM);
ind=collect(1:size(t,1));
ind=ind[i];                                   
#  Take triangles with internal centroids
t = t[ind,:];

nelem=size(ELEM,1);
plt.plot(NOS_GEO[:,1],NOS_GEO[:,2],"kx",markersize=4,linewidth=1);	# Plot the node of the elements
plt.grid("on")

plt.plot(NOS[:,1],NOS[:,2],"ko",markersize=4);	# Plot the node of the elements
ELEM2=[ELEM ELEM[:,2]];

ncont=50;
# plt.triplot(XY[:,1], XY[:,2], t-1, color=(0.0,0.,0.),linewidth=0.4)
plt.triplot(NOS_GEO[:,1], NOS_GEO[:,2], ELEM2-1, color=(0.0,0.,0.),linewidth=0.4)

plt.tricontourf(XY[:,1], XY[:,2], t-1, cor, ncont)
plt.colorbar()

plt.quiver(PONTOS_INT[:,1],PONTOS_INT[:,2],dTdx,dTdy,color="red",width=0.002,scale=50)
plt.axis("equal")
ax = plt.gca() # get current axes
ax[:set_xlim]((xmin-.15*deltax,xmax+.15*deltax));
ax[:set_ylim]((ymin-.15*deltay,ymax+.15*deltay));

end
