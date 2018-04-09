% Entrada de dados para análise de temperatura pelo
% método dos elementos de contorno

% Matriz para definição de pontos que definem a geometria
% PONTO = [número do ponto, coord. x do ponto, coord. y do ponto]

a=1;
b=2;

PONTOS  = [1 -b 0 ;
    2 0 -b ;
    3 b 0 ;
    4 0 b ;
    5 -a 0 ;
    6 0 a ;
    7 a 0
    8 0 -a];

% Segmentos que definem a geometria
%  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
%                                                  Raio]
% Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
%                          inicial para o ponto final)
%                   < 0 -> O centro é à direita do segmento (do ponto
%                          inicial para o ponto final)
%                   = 0 -> O segmento é uma linha reta
SEGMENTOS = [1 1 2 b;
    2 2 3 b;
    3 3 4 b;
    4 4 1 b
    5 5 6 -a
    6 6 7 -a
    7 7 8 -a
    8 8 5 -a];

% Matriz para definição da malha

% MALHA =[numero do segmento, numero de elementos no segmento]
nelem=1;
MALHA = [1 nelem;
    2 nelem;
    3 nelem;
    4 nelem;
    5 nelem;
    6 nelem;
    7 nelem;
    8 nelem];
% Condições de contorno nos segmentos
% CCSeg=[no do segmento, tipo da CDC, valor da CDC]
% tipo da CDC = 0 => a temperatura é conhecida
% tipo da CDC = 1 => O fluxo é conhecido
CCSeg=[1 1 -200
    2 1 -200
    3 1 -200
    4 1 -200
    5 0 100
    6 0 100
    7 0 100
    8 0 100];

% fc = fonte de calor concentrada
% fc = [valor da fonte, coordenada x da fonte, coordenada y da fonte];
fc=[0 .5 .5]; 


k=1;

% Definição dos pontos internos
npi=5;
NPX=npi; % Número de pontos internos na direção X
NPY=npi; % Número de pontos internos na direção Y

