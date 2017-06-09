PONTOS  = [1   0  0
          2    1  0];

% Segmentos que definem a geometria
%  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final
%                                                  Raio, tipo do elemento]
% Raio do segmento: > 0 -> O centro ? ? esquerda do segmento (do ponto
%                          inicial para o ponto final)
%                   < 0 -> O centro ? ? direita do segmento (do ponto
%                          inicial para o ponto final)
%                   = 0 -> O segmento ? uma linha reta


SEGMENTOS = [1 1 2  .5
         2 2 1  .5];

% Matriz para defini??o da malha

% MALHA =[numero do segmento, numero de elementos no segmento]

ne=5;
MALHA = [1 ne
    2 ne];
CCSeg=[1 1 0
    2 0 1];

 kmat=1;
