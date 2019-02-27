function interp(Tree,block,arg,ninterp=3,compressão=true,ϵ=1e-3)
    NOS, NOS_GEO, ELEM, k, CDC,qsi,w,qsi_tri,w_tri = arg;
    #    NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k = arg
    n = size(block,1)               # Number of blocks
    Aaca = Array{Any}(n,2)          
    b = complex(zeros(size(arg[1],1)))       
    for i=1:n                       # Loop over each block
	b1 = Tree[block[i,1]]       
	b2 = Tree[block[i,2]]       
	if block[i,3]==0                # If the block is not admissible for approximation, calculate all the terms
            Aaca[i,1],B = cal_Aeb(b1,b2,arg)
	    b[b1] = b[b1] + B  
	else                            # The block is admissible for approximation
            Aaca[i,1],Aaca[i,2],L,B=cal_Aeb_interp(b1,b2,arg,ninterp,compressão,ϵ)
            b[b1] = b[b1] + L*(B*arg[7][b2,3])
            # Aaca[i,1],B = cal_Aeb(b1,b2,arg)
	    # b[b1] = b[b1] + B  
	end	
    end
    return Aaca,b
end

function cal_Aeb_interp(b1,b2,arg,ninterp=3,compressão=true,ϵ=1e-3)
    NOS, NOS_GEO, ELEM, k, CDC,qsi,w,qsi_tri,w_tri = arg;
    #    NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k=arg
    nelem::Int64 = size(ELEM)[1]          
    G = complex(zeros(ninterp*ninterp*ninterp,length(b2)))      
    H = complex(zeros(ninterp*ninterp*ninterp,length(b2)))      
    xmax=zeros(1,3); xmin=zeros(1,3)
    xmax[1]=maximum(NOS[b1,2])
    xmin[1]=minimum(NOS[b1,2])
    xmax[2]=maximum(NOS[b1,3])
    xmin[2]=minimum(NOS[b1,3])
    xmax[3]=maximum(NOS[b1,4])
    xmin[3]=minimum(NOS[b1,4])
    xs=criapontosinterp(ninterp)
    xint=zeros(ninterp*ninterp*ninterp,3)
    n1,n2=calc_fforma(xs)
    xks=n1 .*xmax .+ n2 .*xmin
    ci=0
    for i3 =1:ninterp # Laco sobre os pontos fontes
        for i2 =1:ninterp # Laco sobre os pontos fontes
            for i1 =1:ninterp # Laco sobre os pontos fontes
                ci+=1
                xd=xks[i1,1]; # Coordenada x do ponto fonte
                yd=xks[i2,2]; # Coordenada y do ponto fonte
                zd=xks[i3,3]; # Coordenada y do ponto fonte
                xyz[ci,:]=[xd yd zd]
                xint[ci,:]=[xd yd zd]
                cj=0       
                for j in b2 # Laco sobre os elementos
                    cj+=1
                    no1::Int64=ELEM[j,1]; # Ponto inicial do elemento
                    no2::Int64=ELEM[j,2]; # Ponto final do elemento
                    no3::Int64=ELEM[j,3]; # Ponto final do elemento
                    x1=NOS_GEO[no1,1]; # Coordenada x do ponto inicial do elemento
                    x2=NOS_GEO[no2,1]; # Coordenada x do ponto final do elemento
                    x3=NOS_GEO[no3,1]; # Coordenada x do ponto final do elemento
                    y1=NOS_GEO[no1,2]; # Coordenada y do ponto inicial do elemento
                    y2=NOS_GEO[no2,2];  # Coordenada y do ponto final do elemento
                    y3=NOS_GEO[no3,2];  # Coordenada y do ponto final do elemento
                    z1=NOS_GEO[no1,3]; # Coordenada z do ponto inicial do elemento
                    z2=NOS_GEO[no2,3];  # Coordenada  do ponto final do elemento
                    z3=NOS_GEO[no3,3];  # Coordenada z do ponto final do elemento
                    n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); 
                    g=0
                    h=0
                    for l =1:npg
                        for m =1:npg
                            qsi = (1 - xi[l])*xi[m]
                            N =calc_fformatri(qsi,xi[l]); #  fun��es de forma
                            x=N[1]*x1+N[2]*x2+N[3]*x3; # coordenada x do ponto de integra��o
                            y=N[1]*y1+N[2]*y2+N[3]*y3 # coordenada y do ponto de integra��o
                            z=N[1]*z1+N[2]*z2+N[3]*z3; # coordenada z do ponto de integra��o
                            J = calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,qsi,xi[m]);# jacobiano
                            Tast,qast=calc_solfund_POT(x,y,z,xd,yd,zd,n,k); # Solu��es
                            g=g+Tast*complex((1-xi[l])*w[l]*w[m]*J,0); # Integral da matriz G
                            h=h+qast*complex((1-xi[l])*w[l]*w[m]*J,0); # Integral da matriz H
                        end
                    end
                    GG[ci,cj] = g
                    HH[ci,cj] = h
                end
            end
        end
    end
    if abs(xmax[1]-xmin[1])<ϵ && abs(xmax[2]-xmin[2])>=ϵ && abs(xmax[3]-xmin[3])>=ϵ
        fontes=[(2. .*(PONTOS_dom[b1,2] .-xmin[2])./(xmax[2]-xmin[2]).-1) (2. .*(PONTOS_dom[b1,3] .-xmin[3])./(xmax[3]-xmin[3]).-1)];
        LL=lagrange(fontes,xs,ninterp,xs,ninterp);
        L=[LL LL]./2;
    elseif abs(xmax[2]-xmin[2])<ϵ && abs(xmax[1]-xmin[1])>=ϵ && abs(xmax[3]-xmin[3])>=ϵ
        fontes=(2. .*(PONTOS_dom[b1,1] .-xmin[1])./(xmax[1]-xmin[1]).-1);
        LL=lagrange(fontes,xs,ninterp,xs,ninterp);
        L=[LL LL]./2;
    elseif abs(xmax[3]-xmin[3])<ϵ && abs(xmax[2]-xmin[2])>=ϵ && abs(xmax[1]-xmin[1])>=ϵ
        fontes=(2. .*(PONTOS_dom[b1,1] .-xmin[1])./(xmax[1]-xmin[1]).-1);
        LL=lagrange(fontes,xs,ninterp,xs,ninterp);
        L=[LL LL]./2;
    elseif abs(xmax[2]-xmin[2])<ϵ && abs(xmax[3]-xmin[3]) < ϵ && abs(xmax[1]-xmin[1])>=ϵ
        fontes=(2. .*(PONTOS_dom[b1,1] .-xmin[1])./(xmax[1]-xmin[1]).-1);
        LL=lagrange(fontes,xs,ninterp);
        L=[LL LL LL]./3;
    elseif abs(xmax[2]-xmin[2])<ϵ && abs(xmax[1]-xmin[1]) < ϵ && abs(xmax[2]-xmin[2])>=ϵ
        fontes=(2. .*(PONTOS_dom[b1,3] .-xmin[3])./(xmax[3]-xmin[3]).-1);
        LL=lagrange(fontes,xs,ninterp);
        L=[LL LL LL]./3;
    elseif abs(xmax[1]-xmin[1])<ϵ && abs(xmax[3]-xmin[3]) < ϵ && abs(xmax[2]-xmin[2])>=ϵ
        fontes=(2. .*(PONTOS_dom[b1,2] .-xmin[2])./(xmax[2]-xmin[2]).-1);
        LL=lagrange(fontes,xs,ninterp);
        L=[LL LL LL]./3;
    else
        fontes=[(2. .*(PONTOS_dom[b1,1] .- xmin[1]) ./(xmax[1]-xmin[1]).-1) (2. .*(PONTOS_dom[b1,2] .-xmin[2])./(xmax[2]-xmin[2]).-1)  (2. .*(PONTOS_dom[b1,3] .-xmin[3])./(xmax[3]-xmin[3]).-1)];
        L=lagrange(fontes,xs,ninterp,xs,ninterp,xs,ninterp);
    end;
L=reverse(L,dims=1);
return L,H,L,G
end
end


#"pg ponto interpolado
#x ponto interpolador"

function lagrange(pg,x,n)
    ni = length(pg);
    L = ones(ni,n);
    for j = 1:n
	for i = 1:n
	    if (i != j)
		L[:,j] = L[:,j].*(pg - x[i])/(x[j]-x[i]);
	    end
	end
    end
    return L
end

function lagrange(pg,x1,n1,x2,n2)
    l1=lagrange(pg[:,1],x1,n1)
    l2=lagrange(pg[:,2],x2,n2)


    ni=size(pg,1)
    L=zeros(ni,n1*n2)
    for i=1:ni
	L[i,:]=(l1[i,:]*l2[i,:]')[:]
    end
    L
end
function lagrange(pg,x1,n1,x2,n2,x3,n3)
    l1=lagrange(pg[:,1:2],x1,n1,x2,n2)
    l2=lagrange(pg[:,3],x3,n3)


    ni=size(pg,1)
    L=zeros(ni,n1*n2*n3)
    for i=1:ni
	L[i,:]=(l1[i,:]*l2[i,:]')[:]
    end
    L
end

function criapontosinterp(n)
    x = cos.((2*(1:n)-1)*pi/2/n);
end
