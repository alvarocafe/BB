% Entrada de dados para análise de temperatura pelo
% método dos elementos de contorno 


% Matriz para definição de pontos que definem a geometria
% PONTO = [número do ponto, coord. x do ponto, coord. y do ponto]

PONTOS  = [1 0 1 
          2 0 0 ;
          3 1 0 ;
          4 1 1 ];

% Segmentos que definem a geometria
%  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
%                                                  Raio] 
% Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
%                          inicial para o ponto final) 
%                   < 0 -> O centro é à direita do segmento (do ponto
%                          inicial para o ponto final)
%                   = 0 -> O segmento é uma linha reta

SEGMENTOS = [1 1 2 0;
	 	   2 2 3 0;
         3 3 4 0;
         4 4 1 0];

% Matriz para definição da malha

% MALHA =[numero do segmento, numero de elementos no segmento]

ne=2;
MALHA = [1  ne;
	       2 ne;
          3 ne;
          4 ne];

   % Condições de contorno nos segmentos
  % CCSeg=[no do segmento, tipo da CDC, valor da CDC]
  % tipo da CDC = 0 => a temperatura é conhecida
  % tipo da CDC = 1 => O fluxo é conhecido
 CCSeg=[1 0 0
     2 1 0
    3 0 1
    4 1 0];

% Condutividade Térmica do material

k = 1;		% [W/m.K]

% fc = fonte de calor concentrada
% fc = [valor da fonte, coordenada x da fonte, coordenada y da fonte];
% fc=[-1 .5 .5];
fc=[0 0 0];

% Definição dos pontos internos
NPX=5;
NPY=5;