function aplica_CDC(G,H,CDC)
    # Applies the boundary conditions by exchanging the columns of matrices G and H

    n_el=size(CDC,1); # Number of boundary conditions
    # Now two new matrices will be built to apply the boundary conditions. The reason for this style of boundary condition collocation is that the boundary condition matrix stored is for the continuous linear element.
    todos_valores=zeros(1,2*n_el);	# Allocate matrix for all the values of the boundary conditions
    todos_tipos = zeros(1,2*n_el);	# Allocate matrix for all the types of boundary conditions
    for i=1:n_el	# Loop over the elements
        todos_valores[2*i-1] = CDC[i,3];
        todos_valores[2*i] = CDC[i,5];
        todos_tipos[2*i-1] = CDC[i,2];
        todos_tipos[2*i] = CDC[i,4];
    end

    ncdc = 2*n_el; # The number of boundary conditions is twice the number of elements, as there are two nodes for each element
    A=copy(H);
    B=copy(G);
    for i=1:ncdc # Loop over the boundary conditions
        tipoCDC = todos_tipos[i]; # Type of boundary condition
        if tipoCDC == 0 # The velocity potential is known
            colunaA=-A[:,i];
            A[:,i]=-B[:,i]; # Matrix A receives the column from matrix G
            B[:,i]=colunaA; # Matrix B receives the column from matrix H
        end
    end;

    b=B*todos_valores'; # vector b
    return A,b
end
