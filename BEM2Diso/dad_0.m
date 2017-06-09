% Entrada de dados 


% Matriz para defini��o de pontos
% PONTO = [n�mero do ponto, coordenada x do ponto, coordenada y do ponto];
L=2;
PONTOS  = [1   0   0 
          2   2*L   0 
          3   2*L  L
			 4   0  L];
          
% Segmentos que definem a geometria
%  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
%                                                  Radio] 
% Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
%                          inicial para o ponto final) 
%                   < 0 -> O centro � � direita do segmento (do ponto
%                          inicial para o ponto final)
%                   = 0 -> O segmento � uma linha reta

SEGMENTOS = [1 1 2 0;
	 	   2 2 3 0;
           3 3 4 0;
	       4 4 1 0];

% Matriz para defini��o da malha
% Malha = [n�mero do segmento, n�mero de elementos no segmento]
ne=4;
MALHA = [1 ne;
	     2 ne;
         3 ne;
		 4 ne];
% Condi��es de contorno nos segmentos
% CCSeg=[Segmento,tipo da CDC em x, valor da CDC em x , ...
%                            tipo da CDC em y, valor da CDC em y]
% tipo da CDC = 0 => o deslocamento � conhecido
% tipo da CDC = 1 => a for�a de superf�cie � conhecida
% Para condi��es de contorno de for�a normal conhecida proceder:
  % tipo da CDC em x = 2, valor da CDC em x = valor da for�a normal
  % Neste caso pode-se atribuir quaisquer valores para tipo da CDC em y e
  % para valor da CDC em y

CCSeg=[1 1 0
    2 0 1
    3 1 0
    4 0 0];

 kmat=1;
