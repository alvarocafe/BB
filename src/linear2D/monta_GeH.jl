function monta_GeH(ELEM,NOS,CDC,k,fc,qsi,w)
# Evaluates matrices H, G and array q, which corresponds to the influences of the velocity potentials, flux and external sources. The discretized boundary integral equation is, in its matricial form, H phi = G qphi + q. 

n_el = length(ELEM[:,1]);	# Number of elements
n_nos = length(NOS[:,1]);	# Number of physical nodes

# Allocating matrices H, G and array q
 H = complex(zeros(n_nos,n_nos));
 G = complex(zeros(n_nos,n_nos));
 q=complex(zeros(n_nos,1));
#H = zeros(n_nos,n_nos)
#G = zeros(n_nos,n_nos)
#q = zeros(n_nos)
# A = complex(zeros(n_nos,n_nos));
# B = complex(zeros(n_nos,n_nos));
for i = 1 : n_nos	# Loop over the source points (physical nodes)
    # Coordinates of the source point
    xd = NOS[i,2];	# x coordinate of the source point
    yd = NOS[i,3];	# y coordinate of the source point
    for j = 1 : n_el	# Loop over the elements
        # Numbering of the nodes of element j
        no1 = ELEM[j,2];	
        no2 = ELEM[j,3];	 
        # Coordinate of the nodes [x1,y1,x2,y2]
        x1 = NOS[no1,2];	y1 = NOS[no1,3];
        x2 = NOS[no2,2];	y2 = NOS[no2,3];

        if ((i == no1) || (i == no2))  # Evaluates the singular integration term when node i belongs to element j
            if(i==no1)	# In this case, the node is at the point described by qsi = -1/2, at a quarter of the way into the line
                eta,Jt = telles(qsi,-1/2);
                g,h =calc_gh_nsing(x1,y1,x2,y2,xd,yd,k,eta,w.*Jt); # Evaluates the singular terms using the kernel for non singular integration
                #h =h;
            else	# In this case, the node is at the point described by qsi = 1/2, at a quarter of the way to the end of the line
                eta,Jt = telles(qsi,1/2);
                g,h =calc_gh_nsing(x1,y1,x2,y2,xd,yd,k,eta,w.*Jt); # Evaluates the singular terms using the kernel for non singular integration
                #h =h;
            end
	else
           # Non singular integration
           g,h = calc_gh_nsing(x1,y1,x2,y2,xd,yd,k,qsi,w); # Evaluates the non singular term
        end
	# Allocates the terms of matrices H and G
        H[i,(2*j-1):2*j] = h;
        G[i,(2*j-1):2*j] = g;
# Applies the boundary conditions
#        if CDC[j,2]==0
#            A[i,2*j-1:2*j] = -h;
#            B[i,2*j-1:2*j] = -g;
#        else
#            A[i,2*j-1:2*j] = g;
#            B[i,2*j-1:2*j] = h;
#        end
    end
    if fc[1,1]!=0
        q[i]=calc_q(xd,yd,fc,k);
    else
        q[i]=0;
    end
end

# Evaluation of the singular terms of matrix H using the rigid body consideration.
for m = 1 : n_nos # Loop over the nodes
    H[m,m] = 0; # Let the diagonal term be zero
    for n = 1 : n_nos	# Loop over the elements
        if n != m	# If the term is not on the diagonal
            H[m,m] = H[m,m] - H[m,n];	# Subtract the term from the diagonal
        end;
    end;
end;

return G,H,q
end
