% Entrada de dados  

% Matriz para definição de pontos 

PONTOS  = [1   0  0 
    2 .5 -.5
          3   1  0
          4 .5 .5];
          
% Segmentos que definem a geometria
%  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
%                                                  Raio] 
% Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
%                          inicial para o ponto final) 
%                   < 0 -> O centro é à direita do segmento (do ponto
%                          inicial para o ponto final)
%                   = 0 -> O segmento é uma linha reta

SEGMENTOS = [1 1 2  .5
     2 2 3  .5
      3 3 4  .5
       4 4 1  .5];

% Matriz para definição da malha

% MALHA =[numero do segmento, numero de elementos no segmento]

MALHA = [1  6
          2 6
          3 6
          4 6];
      
 k=1; % Condutividade térmica do material
 
   % Condições de contorno nos segmentos
  % CCSeg=[no do segmento, tipo da CDC, valor da CDC]
  % tipo da CDC = 0 => a temperatura é conhecida
  % tipo da CDC = 1 => O fluxo é conhecido
 CCSeg=[1 0 0
     2 1 0
     3 1 0
     4 1 0];
 
% Definição dos pontos internos 
NPX=9; % Número de pontos internos na direção X
NPY=9; % Número de pontos internos na direção Y

% fc = fonte de calor concentrada
% fc = [valor da fonte, coordenada x da fonte, coordenada y da fonte];
fc=[-1 .52 .52]; 
