function divnode(X,t)
  n=length(t);
  if n>1
    x=X[t,:];
    c=mean(x,1);
    C=cov(x);
    eig_val,eig_vec=eig(C);
    ind=indmax(eig_val);
    v=real(eig_vec[:,ind]);
    teste=zeros(n);
    for i=1:n
      teste[i]=dot( vec(x[i,:]'-c),v);
    end
    x1=t[teste.>=0];
    x2=t[teste.<0];
    # diam1=2*maximum(teste)/norm(v);
    diam=norm(maximum(x,1)-minimum(x,1))
    # println(diam1-diam)
  else
    x1=0
    x2=0
    diam=0
    c=X[t,:]
  end
  return x1,x2,diam,c
end


function cluster(X, max_elem,η = 1.0)
  #uvec t1,t2,child1(2*max_clt);
  #X=rand(10,2)
  #max_elem=1
  m,n = size(X);
  max_clt = m*10;
  child1 = zeros(1,2*max_clt);
  child2 = zeros(1,2*max_clt);
  t =collect(1:m)
  inode = 1;
  ileaf = 1;
  nodes = Array{Any}(2*max_clt)
  leaves= Array{Any}(2*max_clt)
  child = zeros(Int,2*max_clt,2);
  nodes[1] = t;
  center = zeros(2*max_clt,2);
  #for i=1:inode+1
  center_row=zeros(2*max_clt,2)
  diam=zeros(2*max_clt)
  i = 1;
  while i<inode+1 && inode >= 1
    #[x1,x2,diam,c] = divnode(X,t)
    t1,t2,d,c = divnode(X,nodes[i]);
    center_row[i,:] = c;
    diam[i] = d;
    if length(t1)> max_elem  # is a node and will be further divided
      inode = inode + 1;
      nodes[inode] = t1;
      child[i,1] = inode;
    else # is a leaf, and is the last node of the tree
      leaves[ileaf] = t1;
      ileaf = ileaf + 1;
      child1[i] = ileaf;

    end
    if length(t2) > max_elem
      inode = inode + 1;
      nodes[inode] = t2;
      child[i,2] = inode;
    else
      leaves[ileaf] = t2;
      ileaf = ileaf + 1;
      child2[i] = ileaf;
    end
    i = i + 1;
  end
  #for (unsigned i = 0; i <inode+1; ++i) [

  Tree = Array{Any}(inode+ileaf-1)
  for i=1:inode
    Tree[i] = nodes[i]; # first, all nodes go into the tree
    if child1[i] > 0
      child[i,1] = child1[i] + inode -1; # add last child term in each division
    end
    if child2[i] > 0
      child[i,2] = child2[i] + inode - 1;
    end
  end

  for i=1:ileaf-1
    Tree[inode+i] = leaves[i]; # then add all the leaves !
    t1,t2,d,c = divnode(X,leaves[i]);
    center_row[inode+i,:] = c;
    diam[inode+i] = d;
    child[inode+i,:] = [0 0];
  end

  dist = zeros(inode+ileaf-1,inode+ileaf-1);
  # η = 1.0; # admissibility criteria, no need to be higher than 1
  for i=1:inode+ileaf-1
    for j=1:inode+ileaf-1
      dist[i,j] = η*norm(center_row[i]-center_row[j],2) - 0.5*(diam[i]+diam[j])-min(diam[i],diam[j]);
      #distance between clusters
    end
  end

  allow = dist.>=0; # if 1 then apply ACA, if 0 do not use ACA
  block=blocks(Tree,child,allow)
  return  Tree,block
end

function blocks(Tree,child,allow)
  c=1;
  fc1=[2; 2; 3; 3];
  fc2=[2; 3; 2; 3];
  block=zeros(Int,0,3)
  #for c1=1:length[fc1]/2
  c1 = 0;
  while c1 < length(fc1)/2
    for i=1:2
      if allow[fc1[c1*2+i],fc2[c1*2+i]]==1
        #block[c,1]=fc1[c1*2+i];
        #block[c,2]=fc2[c1*2+i];
        #block[c,3]=1;
        block=vcat(block,[fc1[c1*2+i] fc2[c1*2+i] 1])
        c=c+1;
      else
        if child[fc1[c1*2+i],1]==0 && child[fc2[c1*2+i],1]==0
          #block[c,1]=fc1[c1*2+i];
          #block[c,2]=fc2[c1*2+i];
          #block[c,3]=0;
          block=vcat(block,[fc1[c1*2+i] fc2[c1*2+i] 0])
          c=c+1;
        else
          if length(Tree[fc1[c1*2+i]])>=length(Tree[fc2[c1*2+i]])
            fc1=[fc1; child[fc1[c1*2+i],:]];
            fc2=[fc2; fc2[c1*2+i]; fc2[c1*2+i]];
          else
            fc2=[fc2; child[fc2[c1*2+i],:]];
            fc1=[fc1; fc1[c1*2+i]; fc1[c1*2+i]];
          end
        end
      end
    end
    c1 = c1 + 1;
  end
  #  [fc1;fc2]
  return block
end

function ACAF(Tree,block,fHeG,arg,erro=1e-5)
  n=size(block,1)
  Aaca=Array{Any}(n,2)
  #Baca=Array{Any}(n,2)
  b=complex(zeros(size(arg[1],1)))
  for i=1:n
    b1=Tree[block[i,1]]
    b2=Tree[block[i,2]]
    if block[i,3]==0
      # if 0==0
      Aaca[i,1],B=fHeG(b1,b2,arg)
      b[b1]+=B*arg[9][b2,3]
    else
      # println(i)

      INDB1=[]
      INDB2=[]
      B1=zeros(0,length(b2))
      B2=zeros(length(b1),0)

      ind1=trues(length(b1))
      ind2=trues(length(b2))

      #linha 1
      indaref=1
      aref=0*ind1
      for indaref =1 :length(ind1)
        aref,btemp=fHeG(b1[indaref],b2,arg)
        aref=aref.'
        push!(INDB1,indaref)
        B1=[B1;btemp]
        if norm(aref)>1e-10
          break
        end
        ind1[indaref]=0
      end
      if ind1==falses(ind1)
        Aaca[i,1]=ind1*0
        Aaca[i,2]=ind2.'*0
      else
        #coluna 1
        indbref=1
        bref=0*ind2

        for ii =1 :length(ind2)
          indbref=indmin(abs(aref[ind2]))
          indbref+=cumsum(ind2.==0)[ind2][indbref]
          bref,btemp=fHeG(b1,b2[indbref],arg)
          push!(INDB2,indbref)
          B2=[B2 btemp]
          if norm(bref)>1e-10
            break
          end
          ind2[indbref]=0
        end

        arefmax=indmax(abs(aref))
        brefmax=indmax(abs(bref))
        Umax=0
        Vmax=0
        nmin=min(length(ind1),length(ind2))
        U=complex(zeros(length(ind1),nmin))
        V=complex(zeros(nmin,length(ind2)))
        Ap=complex(zeros(length(ind1),length(ind2)))
        norma0=0.0
        cont=0
        for cont=1:nmin-1
          if abs(aref[arefmax])>abs(bref[brefmax])
            U[:,cont],btemp=fHeG(b1,b2[arefmax],arg)
            U[:,cont]=U[:,cont]-Ap[:,arefmax]
            push!(INDB2,arefmax)
            B2=[B2 btemp]
            Umax=indmax(abs(U[:,cont]))
            V[cont,:],btemp=fHeG(b1[Umax],b2,arg)
            V[cont,:]=(V[cont,:][:]-Ap[Umax,:]).'/U[Umax,cont]
            push!(INDB1,Umax)
            B1=[B1;btemp]
            # V[cont,:]=(fHeG(b1[Umax],b2,arg)[1,:]-Ap[Umax,:]).'/U[Umax,cont]
            Vmax=arefmax
          else
            V[cont,:],btemp=fHeG(b1[brefmax],b2,arg)
            V[cont,:]=(V[cont,:][:]-Ap[brefmax,:])
            push!(INDB1,brefmax)
            B1=[B1;btemp]
            Vmax=indmax(abs(V[cont,:]))
            U[:,cont],btemp=fHeG(b1,b2[Vmax],arg)
            U[:,cont]=(U[:,cont]-Ap[:,Vmax])/V[cont,Vmax]
            push!(INDB2,Vmax)
            B2=[B2 btemp]
            Umax=brefmax
          end
          Ap=U*V
          norma1=vecnorm(Ap)
          if abs((norma1-norma0)/norma1) < erro
            break
          else
            norma0=norma1
          end
          ind1[Umax]=0
          ind2[Vmax]=0
          if indaref==Umax && indbref==Vmax

            for indaref =1 :sum(ind1)
              indaref=indmax(ind1)
              aref,btemp=fHeG(b1[indaref],b2,arg)
#              println("Hello World")
#              println("size(aref)",size(aref))
#              println("size(Ap[indaref,:])",size(Ap[indaref,:]))
              aref=aref-Ap[indaref,:]'
              push!(INDB1,indaref)
              B1=[B1;btemp]
              if norm(aref)>1e-10
                break
              end
              ind1[indaref]=0
            end



            for indbref =1 :sum(ind2)
              indbref=indmin(abs(aref[ind2]))
              indbref=indbref+cumsum(ind2.==0)[ind2][indbref]
              bref,btemp=fHeG(b1,b2[indbref],arg)
              bref=bref-Ap[:,indbref]

              push!(INDB2,indbref)
              B2=[B2 btemp]
              if norm(bref)>1e-10
                break
              end
              ind2[indbref]=0
            end

          elseif indaref==Umax
            bref=bref-U[:,cont]*V[cont,indbref].'
            #  indaref=indmin(abs(bref[ind1]))
            #  indaref=indaref+cumsum(ind1.==0)[ind1][indaref]
            #  aref=fHeG(b1[indaref],b2).'-Ap[indaref,:].'

            for indaref =1 :sum(ind1)
              indaref=indmin(abs(bref[ind1]))
              indaref=indaref+cumsum(ind1.==0)[ind1][indaref]
              aref,btemp=fHeG(b1[indaref],b2,arg)
              aref=aref[:]-Ap[indaref,:][:]
              push!(INDB1,indaref)
              B1=[B1;btemp]
              if norm(aref)>1e-10
                break
              end
              ind1[indaref]=0
            end


          elseif indbref==Vmax
            aref=aref[:]-U[indaref,cont]*V[cont,:]
            # indbref=indmin(abs(aref[ind2]))
            # indbref=indbref+cumsum(ind2.==0)[ind2][indbref]
            # bref=fHeG(b1,b2[indbref])-Ap[:,indbref]

            for indbref =1 :sum(ind2)
              indbref=indmin(abs(aref[ind2]))
              indbref=indbref+cumsum(ind2.==0)[ind2][indbref]
              bref,btemp=fHeG(b1,b2[indbref],arg)
              bref=bref-Ap[:,indbref]

              push!(INDB2,indbref)
              B2=[B2 btemp]
              if norm(bref)>1e-10
                break
              end
              ind2[indbref]=0
            end
          else
            aref=aref[:]-U[indaref,cont]*V[cont,:]
            bref=bref-U[:,cont]*V[cont,indbref].'
          end
          arefmax=indmax(abs(aref))
          brefmax=indmax(abs(bref))
        end
        Aaca[i,1]=U[:,1:cont]
        Aaca[i,2]=V[1:cont,:]
        # Aaca[i,1]=U
        # Aaca[i,2]=V
      end
      max1=ind2sub(size(B1),indmax(abs(B1[:,INDB2])))
      maxv=B1[max1[1],INDB2[max1[2]]]
      Vb=B1[max1[1],:].'
      Ub=B2[:,max1[2]]/maxv
      Bp=Ub*Vb
      B1-=Bp[INDB1,:]
      B2-=Bp[:,INDB2]
      for i=1:size(B1,1)-1
        max1=ind2sub(size(B1),indmax(abs(B1[:,INDB2])))
        maxv=B1[max1[1],INDB2[max1[2]]]
        if abs(maxv)<1e-12
          break
        end
        Vb=[Vb; B1[max1[1],:].']
        Ub=[Ub B2[:,max1[2]]/maxv]
        lastBp=Ub[:,end]*Vb[end,:].'
        B1-=lastBp[INDB1,:]
        B2-=lastBp[:,INDB2]

      end
      b[b1]+=Ub*Vb*arg[9][b2,3]

      # println(Bp)

    end

    # println(INDB2)

  end

  return Aaca,b
end
function erroblocos(hmat,A,block,Tree)
  for i =1:length(block[:,3])
    if block[i,3]==1
      n=vecnorm(A[Tree[block[i,1]],Tree[block[i,2]]]-hmat[i,1]*hmat[i,2])
      println("erro no bloco $i = $n")
    end
  end
end

function matvec(hmat,b,block,Tree)
  v=b*0
  for i =1:length(block[:,3])
    if block[i,3]==1
      v[Tree[block[i,1]]]+=hmat[i,1]*hmat[i,2]*b[Tree[block[i,2]]]
    else
      v[Tree[block[i,1]]]+=hmat[i,1]*b[Tree[block[i,2]]]
    end
  end
  v
end
function montacheia(hmat,block,Tree,n)
A=complex(zeros(n,n))
  for i =1:length(block[:,3])
    if block[i,3]==1
      A[Tree[block[i,1]],Tree[block[i,2]]]=hmat[i,1]*hmat[i,2]
    else
      A[Tree[block[i,1]],Tree[block[i,2]]]=hmat[i,1]
    end
  end
  A
end
