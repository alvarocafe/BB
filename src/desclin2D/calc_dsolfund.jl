function  [dTdx,dTdy,dqdx,dqdy]=calc_dsolfund(x,y,xd,yd,nx,ny,k) 
%Calcula as soluções fundamentais

r=sqrt((x-xd)^2+(y-yd)^2); % Raio (distância entre ponto fonte e 
                           % ponto campo)
rx=(x-xd); % Componente x do raio
ry=(y-yd); % Componente y do raio

dTdx=1/(2*pi*k*r^2)*rx; % Solução fundamental da temperatura
dTdy=1/(2*pi*k*r^2)*ry; % Solução fundamental da temperatura

dqdx=(nx*(rx^2 - ry^2) + 2*ny*rx*ry)/(2*pi*r^4);
dqdy=(ny*(-rx^2 +ry^2) + 2*nx*rx*ry)/(2*pi*r^4);