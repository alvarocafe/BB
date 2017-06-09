function [T,q]= monta_Teq(CDC,x)
% Separa fluxo e temperatura

% ncdc = número de linhas da matriz CDC
% T = vetor que contém as temperaturas nos nós
% q = vetor que contém o fluxo nos nós

ncdc = length(CDC(:,1));

for i=1:ncdc % Laço sobre as condições de contorno
    tipoCDC=CDC(i,2); % Tipo da condição de contorno
    valorCDC=CDC(i,3); % Valor da condição de contorno
    valorcalculado=x(i); % Valor que antes era desconhecido
    if tipoCDC == 1 % Fluxo é conhecido
        T(i) = valorcalculado; % A temperatura é o valor calculado
        q(i) = valorCDC; % O fluxo é a condiçao de contorno
    else % A temperatura é conhecida
        T(i) = valorCDC; % A temperatura é a condiçao de contorno
        q(i) = valorcalculado; % O fluxo é o valor calculado
    end
end

