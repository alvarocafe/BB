function Monta_Teq(x,CDC)

# Function to reorder de vectors in accordance with boundary conditions
# Separa fluxo e temperatura

# ncdc = n�mero de linhas da matriz CDC
# T = vetor que cont�m as temperaturas nos n�s
# q = vetor que cont�m o fluxo nos n�s

n_el=size(CDC,1);
todos_valores=zeros(1,2*n_el);
todos_tipos = zeros(1,2*n_el);
for i=1:n_el
    todos_valores[2*i-1] = CDC[i,3];
    todos_valores[2*i] = CDC[i,5];
    todos_tipos[2*i-1] = CDC[i,2];
    todos_tipos[2*i] = CDC[i,4];
end

ncdc = 2*n_el;

# T = complex(zeros(ncdc))
# q = complex(zeros(ncdc))
T = zeros(ncdc)
q = zeros(ncdc)
for i=1:ncdc # La�o sobre as condi��es de contorno
    tipoCDC=todos_tipos[i]; # Tipo da condi��o de contorno
    valorCDC=todos_valores[i]; # Valor da condi��o de contorno
    valorcalculado=x[i]; # Valor que antes era desconhecido
    if tipoCDC == 1 # Fluxo � conhecido
        T[i] = valorcalculado; # A temperatura � o valor calculado
        q[i] = valorCDC; # O fluxo � a condi�ao de contorno
    else # A temperatura � conhecida
        T[i] = valorCDC; # A temperatura � a condi�ao de contorno
        q[i] = valorcalculado; # O fluxo � o valor calculado
    end
end


return T,q
end
