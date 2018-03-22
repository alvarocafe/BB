% Entrada de dados para análise de temperatura pelo
% método dos elementos de contorno (linear)
%
%   Autor: Frederico Lourenço
%   Data de criação junho de 1999 
%   Revisão 0.0
%

% Matriz para definição de pontos 

PONTOS  = [1    0   0 
          2   10   0 
          3   10  10 
			 4    0  10
          5    7   4
          6    4   4
			 7    4   6
			 8    7.5   6];
          
% Segmentos que definem a geometria
%  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
%                                                  Raio] 
% Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
%                          inicial para o ponto final) 
%                   < 0 -> O centro é à direita do segmento (do ponto
%                          inicial para o ponto final)
%                   = 0 -> O segmento é uma linha reta
SEGMENTOS = [1 1 2 0 
	 	   2 2 3 0
         3 3 4 0
			4 4 1 0
			5 5 6 0
         6 6 7 0
			7 7 8 0
			8 8 5 0];

% Matriz para definição da malha

% MALHA =[numero do segmento, numero de elementos no segmento]

MALHA = [1  8 
	       2  8 
          3  8
			 4  8
			 5  6
			 6  6
			 7  6
   		 8  6];

   % Condições de contorno nos segmentos
  % CCSeg=[no do segmento, tipo da CDC, valor da CDC]
  % tipo da CDC = 0 => a temperatura é conhecida
  % tipo da CDC = 1 => O fluxo é conhecido
 CCSeg=[1 1 0
    2 0 0
    3 1 0
    4 0 50
    5 1 -50
    6 1 -10
    7 1 -10
    8 1 -30];

% Condutividade Térmica do material

k = 1;		% [W/m.K]

% fc = fonte de calor concentrada
% fc = [valor da fonte, coordenada x da fonte, coordenada y da fonte];
fc=[-1 .5 .5]; 

% Definição dos pontos internos
NPX=14;
NPY=14;