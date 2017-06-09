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

kmat=1;

% Geração automática de pontos internos
npi=17;
NPX=npi; % Número de pontos internos na direção X
NPY=npi; % Número de pontos internos na direção Y

% disp('Erros em relação a solução analítica')
% nnos=size(NOS,1);
% npint=size(PONTOS_int,1);
% r=sqrt(PONTOS_int(:,2).^2+PONTOS_int(:,3).^2);
a=1;
b=2;
Ti=100;
qo=-200;
qi=-qo*b/a;
To=Ti-qo*b*log(b/a);
% Ta=Ti-qo*b*log(r/a);
% qa=-qo*b./r;
% qr=sqrt(qx.^2+qy.^2);
% [T_pint Ta]
% [qr qa]
% % ErroRMS=norm(T-Ta)/sqrt(npint)
% % Erromax=norm(T-Ta,inf)
% ErroRMS=norm(T_pint-Ta)/sqrt(npint)
% Erromax=norm(T_pint-Ta,inf)