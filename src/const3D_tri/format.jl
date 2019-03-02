function lermsh(nome,dim=2,CCSeg=false)
    dados=readdlm(nome)
    index=findn(dados.=="\$Nodes")
    if(isempty(index[1]))
      index=findn(dados.=="\$ParametricNodes")
  end
  npontos=dados[index[1][1]+1,1]
  pontos=Array{Float64}(dados[index[1][1]+2:index[1][1]+npontos+1,2:4])
  index[1][1]=index[1][1]+4+npontos
  nelem=dados[index[1][1],1]
  elemtipo=dados[index[1][1]+1:index[1][1]+nelem,2]
  elemint=Array{Int}[]
  if dim==2
    if sum(elemtipo.==8)!= 0 # 3-node second order line
        elemcont=Array{Int}(dados[index[1][1]+find(elemtipo.==8),[6:8;5]])
        # seg=Array{Int}(dados[index[1][1]+find(elemtipo.==8),4])
    elseif sum(elemtipo.==1)!= 0 # 2-node line.
        elemcont=Array{Int}(dados[index[1][1]+find(elemtipo.==1),[6:7;5]])
        # seg=Array{Int}(dados[index[1][1]+find(elemtipo.==1),4])
    end

elseif dim==3
    if sum(elemtipo.==2)!= 0 #3-node triangle.
        elemcont=Array{Int}(dados[index[1][1]+find(elemtipo.==2),[6:8;5]])
        elemnum=Array{Int}(dados[index[1][1]+find(elemtipo.==2),1])
        # seg=Array{Int}(dados[index[1][1]+find(elemtipo.==2),4])
    elseif sum(elemtipo.==9)!= 0 #6-node second order triangle
        elemcont=Array{Int}(dados[index[1][1]+find(elemtipo.==9),[6:11;5]])
        # seg=Array{Int}(dados[index[1][1]+find(elemtipo.==9),4])
    end
    if sum(elemtipo.==4)!= 0 #4-node tetrahedron.
        elemcont=Array{Int}(dados[index[1][1]+find(elemtipo.==4),[6:9;5]])
        elemint=Array{Int}(dados[index[1][1]+find(elemtipo.==4),[6:9;5]])
    elseif sum(elemtipo.==11)!= 0 #10-node second order tetrahedron
        elemcont=Array{Int}(dados[index[1][1]+find(elemtipo.==11),[6:15;5]])
        elemint=Array{Int}(dados[index[1][1]+find(elemtipo.==11),[6:15;5]])
    elseif sum(elemtipo.==3)!= 0 #4-node quadrangle
        elemcont=Array{Int}(dados[index[1][1]+find(elemtipo.==3),[6:9;5]])    	
    end

end
if CCSeg==false
    # return pontos,elemcont,elemint
    return [1:size(pontos,1) pontos],[1:size(elemcont,1) elemcont],[1:size(elemint,1) elemint],elemtipo


else
   CDC=Array{Float64}(elemcont)
   for i=1:size(CCSeg,1)
     index=find(seg.==i)
     CDC[index,:]=CCSeg[i,:]
 end
 return pontos,elemcont,elemint,CDC
end
end

function salva(nome,T,elemtipo,tipo=2)
    dados=readdlm(nome)
    index=findfirst(dados[:,1].=="\$ElementData")
        if isempty(index)
        open(nome, "a") do f
            writedlm(f,["\$ElementData"
                1
                "Temperatura"
                1
                0.0
                3
                0
                1
                size(T,1)])
            writedlm(f,[(1:size(T,1))+findfirst(elemtipo,tipo)-1 T])
            write(f,"\$EndElementData")
        end
        else
            open(nome, "w") do f
            writedlm(f,dados[1:index-1,:])
            writedlm(f,["\$ElementData"
                1
                "Temperatura"
                1
                0.0
                3
                0
                1
                size(T,1)])
            writedlm(f,[string.(collect((1:size(T,1))+findfirst(elemtipo,tipo)-1)) T])
            write(f,"\$EndElementData")
        end
    end
end
function mostra_geo(NOS_GEO,ELEM)
# geometria da estrutura

nelem=size(ELEM,1);
xmin = mini(NOS_GEO[:,2]);
xmax = maxi(NOS_GEO[:,2]);
ymin = mini(NOS_GEO[:,3]);
ymax = maxi(NOS_GEO[:,3]);
zmin = mini(NOS_GEO[:,4]);
zmax = maxi(NOS_GEO[:,4]);
NOS=zeros(nelem,4);

for i = 1:nelem
    no1=ELEM[i,2];
    no2=ELEM[i,3];
    no3=ELEM[i,4];
    no4=ELEM[i,5];

    xno1=NOS_GEO[no1,2];
    yno1=NOS_GEO[no1,3];
    zno1=NOS_GEO[no1,4];

    xno2=NOS_GEO[no2,2];
    yno2=NOS_GEO[no2,3];
    zno2=NOS_GEO[no2,4];

    xno3=NOS_GEO[no3,2];
    yno3=NOS_GEO[no3,3];
    zno3=NOS_GEO[no3,4];

    xno4=NOS_GEO[no4,2];
    yno4=NOS_GEO[no4,3];
    zno4=NOS_GEO[no4,4];

    xnos=[xno1,xno2,xno3,xno4];
    ynos=[yno1,yno2,yno3,yno4];
    znos=[zno1,zno2,zno3,zno4];
    #xlabel ('{\it x}')
    #ylabel ('{\it y}')
    #zlabel ('{\it z}')
    plot3D(xnos,ynos,znos,linestyle="none",marker="o");	#Parte do programa que mostra geometria
    #axis equal
    #axis([xmin xmax ymin ymax zmin zmax])
    #hold on;
    xno=sum(xnos)/4;
    yno=sum(ynos)/4;
    zno=sum(znos)/4;
    #text(xno,yno,zno,num2str(i))
    NOS[i,:]=[i,xno,yno,zno];
end;
#view(3)
#hold off;

return NOS
end

function mostra_geoTRI(NOS_GEO,ELEM)
# geometria da estrutura

nelem=size(ELEM,1);
# xmin = min(NOS_GEO(:,2));
# xmax = max(NOS_GEO(:,2));
# ymin = min(NOS_GEO(:,3));
# ymax = max(NOS_GEO(:,3));
# zmin = min(NOS_GEO(:,4));
# zmax = max(NOS_GEO(:,4));
NOS=zeros(nelem,4);

for i = 1:nelem
		no1::Int=ELEM[i,2];
		no2::Int=ELEM[i,3];
		no3::Int=ELEM[i,4];

		xno1=NOS_GEO[no1,2];
		yno1=NOS_GEO[no1,3];
		zno1=NOS_GEO[no1,4];

		xno2=NOS_GEO[no2,2];
		yno2=NOS_GEO[no2,3];
		zno2=NOS_GEO[no2,4];

		xno3=NOS_GEO[no3,2];
		yno3=NOS_GEO[no3,3];
		zno3=NOS_GEO[no3,4];

		xnos=[xno1,xno2,xno3];
    ynos=[yno1,yno2,yno3];
    znos=[zno1,zno2,zno3];
#     xlabel ('{\it x}')
#     ylabel ('{\it y}')
#     zlabel ('{\it z}')
#     plot3D(xnos,ynos,znos,linestyle="none",marker="o");	#Parte do programa que mostra geometria
#     axis equal
#     axis([xmin xmax ymin ymax zmin zmax])
#     hold on;
    xno=sum(xnos)/3;
    yno=sum(ynos)/3;
    zno=sum(znos)/3;
#     text(xno,yno,zno,num2str(i))
    NOS[i,:]=[i,xno,yno,zno];
end;
# view(3)
# hold off;
return NOS
end

function mostra_resultados2(XYZ,tri,T)
nelem=size(tri,1)
x=zeros(3)
y=zeros(3)
z=zeros(3)
zc=zeros(nelem,1)
pc=[zeros(3,3)];
triang=zeros(3,3)
for elem =1:nelem
    no1=tri[elem,2]
    no2=tri[elem,3]
    no3=tri[elem,4]
    x[1]=XYZ[no1,2]
    y[1]=XYZ[no1,3]
    z[1]=XYZ[no1,4]
    x[2]=XYZ[no2,2]
    y[2]=XYZ[no2,3]
    z[2]=XYZ[no2,4]
    x[3]=XYZ[no3,2]
    y[3]=XYZ[no3,3]
    z[3]=XYZ[no3,4]
    triang=[[x[1] y[1] z[1]
    x[2] y[2] z[2]
    x[3] y[3] z[3]]]
    append!(pc,triang)
end
fig = plt.figure()
ax = mp.Axes3D(fig)
q = ar.Poly3DCollection(pc[2:end], linewidths=1)
ax[:add_collection3d](q)
m = cm.ScalarMappable(cmap=cm.jet)
b=m[:to_rgba](T[1:nelem])
q[:set_facecolor](b[:,1:3])
m[:set_array]([minimum(T),maximum(T)])
m[:set_clim](vmin=minimum(T),vmax=maximum(T))
plt.colorbar(m, orientation="vertical",shrink=0.9)
ax[:set_xlabel]("x")
ax[:set_ylabel]("y")
ax[:set_zlabel]("z")
ax[:set_xlim](minimum(XYZ[:,2]),maximum(XYZ[:,2]))
ax[:set_ylim](minimum(XYZ[:,3]),maximum(XYZ[:,3]))
ax[:set_zlim](minimum(XYZ[:,4]),maximum(XYZ[:,4]))
ax[:set_aspect]("equal")
ax[:view_init](elev=18., azim=43.)
return
end
