function cal_GeH(NOS,ELEM,CW,FR,fc,qsi,w)
# Evaluates the G and H matrices for the linear system H phi = G q, where phi is a vector containing the values of the velocity potential and q is a vector containing the values of the flux at the boundary of the problem.
 
  nelem::Int8=size(ELEM)[1]; # Number of elements
  nnos::Int8=size(NOS)[1]; # Number of nodes
  G=complex(zeros(nnos,nnos)); 	# Allocates matrix G
  H=complex(zeros(nnos,nnos));	# Allocates matrix H
  q=zeros(nnos,1);  # Influence from concentrated sources
  inc = zeros(nnos,1);  # Influence from incident plane waves
  qsitelles1,Jtelles1 = telles(qsi,-1); # Evaluates the Telles' points and Jacobian if the field point corresponds to the first node of the element
  qsitelles2,Jtelles2 = telles(qsi,0); # Evaluates the Telles' points and Jacobian if the field point corresponds to the second node of the element
  qsitelles3,Jtelles3 = telles(qsi,1); # Evaluates the Telles' points and Jacobian if the field point corresponds to the third node of the element

  for i=1:nnos # Loop over the source points
    xd=NOS[i,2]; # x coordinate of the source point
    yd=NOS[i,3]; # y coordinate of the source point
    for j=1:nelem # Loop over the elements
      no1::Int8=ELEM[j,2]; # First point of the element
      no2::Int8=ELEM[j,3]; # Second point of the element
      no3::Int8=ELEM[j,4]; # Second point of the element
      x1=NOS[no1,2]; # x coordinate of the first point of the element
      x2=NOS[no2,2]; # x coordinate of the second point of the element
      x3=NOS[no3,2]; # x coordinate of the second point of the element
      y1=NOS[no1,3]; # y coordinate of the first point of the element
      y2=NOS[no2,3];  # y coordinate of the second point of the element
      y3=NOS[no3,3];  # y coordinate of the second point of the element
      if i==no1 # The source point belongs to the element
        g, h = calcula_GeHns(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsitelles1,w.*Jtelles1,FR);	# Singular integration using the Telles transformation
        #h = h - 0.5	# Adding the jump term
        elseif i==no2 # The source point doesnt belong to the element
	g, h = calcula_GeHns(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsitelles2,w.*Jtelles2,FR);	# Singular integration using the Telles transformation
        #h = h - 0.5	# Adding the jump term
	elseif i==no3
	g, h = calcula_GeHns(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsitelles3,w.*Jtelles3,FR);	# Singular integration using the Telles transformation      
        #h = h - 0.5	# Adding the jump term
	else
	g, h = calcula_GeHns(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsi,w,FR);	#Regular integration using Gaussian quadrature
	end
	println("Tipo e tamanho de g:\n tipo=",typeof(g),",\n tamanho=",size(g),",\n valor g=",g)
	println("Tipo e tamanho de h:\n tipo=",typeof(h),",\n tamanho=",size(h),",\n valor h=",h)
	# allocate evaluated integrals into the global H and G matrices
 	for nolocal=1:3
            noglobal::Int8=ELEM[j,nolocal+1]; #Global index of matrix H
            H[i,noglobal] = H[i,noglobal] + h[nolocal];
        end
	println("Elemento = $j")
            G[i,3*j-2:3*j] = g
    end
  end
  # Rigid body consideration for the evaluation of the singular terms of matrix H
for m = 1 : nnos
    H[m,m] = 0; # let the diagonal be zero
    for n = 1 : nnos	# Loop over the source points
        if n != m	# If it's not a diagonal term
            H[m,m] = H[m,m] - H[m,n]; # Subtracts the term from the diagonal
        end
    end
end

  return G,H
end
function Monta_Teq(CDC,ELEM,x)

nelem=length(ELEM[:,1]); #Número de elementos
tipoCDCultimo=CDC[nelem,2];
n_pontos=length[x];
T = complex(zeros())
for el = 1:nelem  # Corre os elementos
    no1 = ELEM[el,2]; # Número do primeiro nó do elemento
    no2 = ELEM[el,3]; # Número do segundo nó do elemento
    no3 = ELEM[el,4]; # Número do segundo nó do elemento
    tipoCDC=CDC[el,2]; # Tipo da condição de contorno no elemento ij
    valorCDC=CDC[el,3]; # Valor da condição de contorno no elemento ij
    if tipoCDC==0   # Temperatura conhecida
        T[[no1,no2,no3]]=[valorCDC valorCDC valorCDC]; # Atribui os valores
        # das condições de contorno ao vetor T
        q[el,1:4]=[el x[no1] x[no2] x[no3]]; # Atribui os valores calculados
        q_vet[3*el-2:3*el]=[x[no1] x[no2] x[no3]];
        # ao vetor q
    else # fluxo conhecido
        q[el,1:4]=[el valorCDC valorCDC valorCDC]; # Atribui os valores
        # das condições de contorno ao vetor q
        q_vet[3*el-2:3*el]=[valorCDC valorCDC valorCDC];
        if el==1 && tipoCDCultimo==1
            T([no1 no2])=x([no1 no2 no3]);
        elseif(el~=nelem)
            T([no2 no3])=[x[no2] x[no3]];
        else
            T[no2]=x[no2];
        end
    end
end
T[2*nelem+1:n_pontos]=x[2*nelem+1:n_pontos];
return T,q,q_vet
end
function calcula_GeHns(x1,y1,x2,y2,x3,y3,xd,yd,CW,qsi,w,FR)
# Non singular integration
  n_pint=size(qsi,1); # Number of integration points
  h = complex(zeros(3)) # Initializes the matrix H terms
  g = complex(zeros(3)) # Initializes the matrix G terms

  for kk=1:n_pint # Loop over the integration points
  N=calc_fforma(qsi[kk]); # Evaluates the shape functions
  dN= calc_dfforma(qsi[kk]) # Evaluates the derivative of the shape functions
  dgamadqsi=cal_Jacobiano([x1 x2 x3],[y1 y2 y3],dN); # Evaluates the Jacobian
  xx=N[1]*x1+N[2]*x2+N[3]*x3; # x coordinate of the integration point
  yy=N[1]*y1+N[2]*y2+N[3]*y3; # y coordinate of the integration point
  dx=dN[1]*x1+dN[2]*x2+dN[3]*x3; # x coordinate derivative of the curve at the integration point
  dy=dN[1]*y1+dN[2]*y2+dN[3]*y3; # x coordinate derivative of the curve at the integration point
  sx=dx/dgamadqsi # x component of the tangent vector
  sy=dy/dgamadqsi # y component of the tangent vector
  nx=sy; # x component of the normal vector
  ny=-sx; # y component of the normal vector
#  phiast,qast =calc_solfund(xx,yy,xd,yd,nx,ny,CW,FR); # Evaluation of the fundamental solutions
	phiast = 1	#uncomment these lines to obtain the length of the boundary
	qast = 1
  g=g+N'*phiast*w[kk]*dgamadqsi # Evaluation of integral for g term
  h=h+N'*qast*w[kk]*dgamadqsi # Evaluation of integral for h term
  end
  return g,h
end

function aplica_CDC(G,H,CDC,ELEM)
# Aplica as condições de contorno (troca as colunas das matrizes G e H)

# Cria a variável T_PR que contém os nós onde a temperatura é prescrita
# (conhecida).
# T_PR tem 5 colunas e o número de linhas é igual ao número de nós onde a
# temperatura é conhecida.
# T_PR=[a1,a2,a3,a4,a5]
# a1=número do nó com temperatura prescrita.
# a2=número do primeiro elemento com temperatura prescrita ao qual este nó
# pertence.
# a3=número local do nó neste elemento.
# a4: caso a temperatura também seja prescrita no segundo elemento a que
# este nó pertence, então, a4 conterá o número deste elemento, caso
# contrário, conterá zero.
# a5: caso a4 seja diferente de zero, a5 conterá o número local do nó no
# segundo elemento, caso contrário, conterá zero.
n_el=size(CDC,1);
todos_valores=zeros(1,3*n_el);
# Determinar o tamanho de T_PR
iter = 0;
for i = 1:size(CDC,1)
	if CDC[i,2*no] == 0
		iter +=1;
	end
end
T_PR=zeros(iter,5);
cont=0;


for el=1:n_el # for sobre os elementos para criar T_PR
    for no=1:3 # for sobre os nós do elemento el
        no_global=ELEM[el,no+1]; # número global do nó
        tipoCDC=CDC[el,2*no];  # tipo da condição de contorno
        valorCDC=CDC[el,2*no+1]; # valor da condição de contorno
        todos_valores[3*el-3+no]=valorCDC; # armazerna o valor da condição
                   # de contorno no vetor todos_valores
        if(tipoCDC==0) # se a temperatura é conhecida
            compartilha=0; # por enquanto não se sabe se a tempertura é 
            # também conhecida no segundo elemento a que este nó pertence
            i=1;
            while (!compartilha&&i<=cont&&cont>0) # verifica se o nó global
                # já está presente em T_PR. Quando ele encontra,
                # compartilha se torna igual a 1 e o while pára.
                if(no_global==T_PR[i,1]) # Se sim, o nó global já está 
                    # presente em T_PR. As colunas 4 e 5 são preenchidas e
                    # compartilha se torna igual a 1.
                    compartilha=1;
                    T_PR[i,4]=el;
                    T_PR[i,5]=no;
                end
                i=i+1;
            end
            if(!compartilha) # Se compartilha continua zero, então o nó 
                # ainda não foi inserido em T_PR. Neste caso, as três
                # primeiras colunas de T_PR são preenchidas.
                cont=cont+1;
                T_PR[cont,1]=no_global;
                T_PR[cont,2]=el;
                T_PR[cont,3]=no;
            end
        end
    end
end

n_temp_pr = length(T_PR[:,1]); # Número de nós com temperatura conhecidas
for i=1 : n_temp_pr # for sobre os nós com temperatura conhecida
    i_no=T_PR[i,1]; # número do nó com temperatura conhecida
    i_el=T_PR[i,2]; # primeiro elemento que contém o nó
    i_no_loc=T_PR[i,3]; # número local do nó neste elemento
    ind_H=i_no; # índice da coluna da matriz H que será trocada
    ind_G=3*i_el+i_no_loc-3; # índice da coluna da matriz G que será 
                 #  trocada
    troca = G[:,ind_G]; # armazena a coluna de G que será trocada
    G[:,ind_G] = -H[:,ind_H]; # substitui a coluna de H na de G
    H[:,ind_H] = -troca; # Substitui a de G na de H
    if(T_PR[i,4]!=0) # Se a temperatura também é conhecida no segundo 
         #   elemento
        i_el=T_PR[i,4]; # Número do segundo elemento 
        i_no_loc=T_PR[i,5]; # número local deste nó no segundo elemento
        ind_G=3*i_el+i_no_loc-3; # índice da coluna G que será atribuído 
                  # zero
        H[:,ind_H]=H[:,ind_H]-G[:,ind_G]; # soma o valor da coluna G na 
                  # coluna H
        G[:,ind_G]=0; # atribui zero na coluna G
    end;
end;
b=G*todos_valores'; # Cálculo do vetor b
A=H;
return A,b,T_PR
end
function calc_solfundpot(x,y,xd,yd,nx,ny,tmp,k)
# Evaluates the fundamental solutions of the Laplace equation.

r=sqrt((x-xd)^2+(y-yd)^2); # Distance between the source and field points
rx=(x-xd); # x component of the distance
ry=(y-yd); # y component of the distance
Tast=-1/(2*pi*k)*log(r); # Fundamental solution for the temperature
qast=1/(2*pi)*(rx*nx+ry*ny)/r^2; # Fundamental solution for the flux
return Tast, qast
end

function  calc_solfund(x,y,xd,yd,nx,ny,CW,FR)
# Evaluates the fundamental solutions of the Helmholtz equation.
r=sqrt((x-xd)^2+(y-yd)^2); # Distance between the source and field points
rx=(x-xd); # x component of the distance
ry=(y-yd); # y component of the distance
drdn=rx*nx+ry*ny;   # Distance in the normal direction
ZR=real(FR*r/CW);
Z=complex(0.,ZR);
F0C=SpecialFunctions.besselk(0,Z);
F1C=SpecialFunctions.besselk(1,Z);

  qast=-(Z/r*drdn*F1C)/(2*pi); 	# Fundamental solution for the velocity potential
  Tast=F0C/(2*pi);    		# Fundamental solution for the flux
  return Tast,qast
end


