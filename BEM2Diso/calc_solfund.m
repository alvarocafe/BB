function  [Tast,qast]=calc_solfund(xd,yd,x,y,nx,ny,k) 
%Calcula as soluções fundamentais

r=sqrt((x-xd)^2+(y-yd)^2); % Raio (distância entre ponto fonte e 
                           % ponto campo)
rx=(x-xd); % Componente x do raio
ry=(y-yd); % Componente y do raio
Tast=-1/(2*pi*k)*log(r); % Solução fundamental da temperatura
qast=1/(2*pi)*(rx*nx+ry*ny)/r^2; % Solução fundamental do fluxo
% disp([r,rx,ry,nx,ny,Tast,qast])
return