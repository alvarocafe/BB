function mapa_de_cor(temp,dTdx,dTdy,NOS,ELEM,PONTOS_INT,PONTOS,SEGMENTOS,T_in)
% Rotina para criação do mapa de cores para temperatura
%
%   Autor: Frederico Lourenço
%   Data de criação julho de 1999
%   Revisão 0.0

% Atribuição das coordenadas e temperatura nos nós e pontos
% internos à matriz Tf[x,y,T]
n_nos=length(ELEM(:,1));
n_pi=length(PONTOS_INT(:,1));
n_linhas = length(SEGMENTOS(:,1));
xmin = min(PONTOS(:,2));
xmax = max(PONTOS(:,2));
lx = xmax - xmin;		% Largura do retângulo que contém a geometria

Tf = zeros(n_nos+n_pi,3);

% Atribuição das coordenadas e temperatura nos nós e pontos
% internos à matriz Tf[x,y,T]
for n = 1 : n_nos		% Percorre todos os nós
  Tf(n,:) = [NOS(n,2),NOS(n,3), temp(n)];
end;

for n = 1 : n_pi		% Percorre todos os pontos internos
  xpi = PONTOS_INT(n,2);
  ypi = PONTOS_INT(n,3);
  Tf(n_nos + n,:) = [xpi,ypi,T_in(n)];
end;

% Geração de cada triângulo da malha com interpolação bilinear
% dos resultados obtidos
for t = 1 : n_linhas		% Percorre todas as linhas
  xl1(t) = PONTOS(SEGMENTOS(t,2),2);
  xl2(t) = PONTOS(SEGMENTOS(t,3),2);
  yl1(t) = PONTOS(SEGMENTOS(t,2),3);
  yl2(t) = PONTOS(SEGMENTOS(t,3),3);
  raio(t)= SEGMENTOS(t,4);  
end;


% Criação dos vetores dos valores de x, de y e de T
xi = Tf(:,1);
yi = Tf(:,2);
ci = Tf(:,3);

% Chamada da função de triangularização (geração de malha
% na domínio para mostrar o mapa de cores)
tri = delaunay(xi,yi);

% Geração de cada triângulo da malha com interpolação bilinear
% das temperaturas
for t = 1 : length(tri)	% Percorre todos os triângulos da malha

  x = [xi(tri(t,1));xi(tri(t,2));xi(tri(t,3))];
  y = [yi(tri(t,1));yi(tri(t,2));yi(tri(t,3))];
  c = [ci(tri(t,1));ci(tri(t,2));ci(tri(t,3))];
  
  % Verificando se o triângulo pertence ao domínio
  xm = (x(1)+x(2)+x(3))/3;
  ym = (y(1)+y(2)+y(3))/3;
  [xm,ponto] = testa_ponto(xm,ym,xl1,yl1,xl2,yl2,lx,raio);

  if strcmp(ponto,'interno')
    patch(x,y,c)
  end;

end;
hold on
view(0,90)
colorbar 			% Cria a barra de cores
title('Distribuição de temperatura');
shading interp
quiver(PONTOS_INT(:,2),PONTOS_INT(:,3),dTdx,dTdy)
hold off
