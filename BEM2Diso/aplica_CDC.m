function [A,b]= aplica_CDC(G,H,CDC,E)
% Aplica as condições de contorno trocando as colunas das matrizes H e G

ncdc = length(CDC(:,1)); % número de linhas da matriz CDC
A=H;
B=G;
for i=1:ncdc % Laço sobre as condições de contorno
    tipoCDC = CDC(i,2); % Tipo da condição de contorno
    if tipoCDC == 0 % A temperatura é conhecida
        colunaA=-A(:,i); % Coluna da matriz H que será trocada
        A(:,i)=-B(:,i); % A matriz H recebe a coluna da matriz G
        B(:,i)=colunaA; % A mstriz G recebe a coluna da matriz H
    end
end

valoresconhecidos=E\CDC(:,3); % Valores das condições de contorno
b=B*valoresconhecidos; % vetor b
