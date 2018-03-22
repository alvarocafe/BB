function [xpi,ponto] = testa_ponto(xpi,ypi,xl1,yl1,xl2,yl2,lx,raio)
% Rotina para verificação da condição de um ponto com relação
% à geometria (interno ou externo).
% O algoritmo utilizado é baseado no número de intersecções entre
% o contorno e uma linha vertical que parte do ponto interno.
% N. ímpar de intersecções - ponto interno
% N. par de intersecções   - ponto externo
%
%   Autor: Frederico Lourenço
%   Data de criação setembro de 1999 
%   Revisão 0.0


pc = 0;					% Contador de pontos do contorno
npc = length(xl1);	% N. de pontos do contorno
sai = 'não';
		
while (strcmp(sai,'não') && pc < npc)
  pc = pc + 1;
  if xpi == xl1(pc)
    xpi = xpi + lx*10^(-2);
    sai = 'sim';
  end;
end;

interv = 0;		% N. de intersecções entre o contorno e a linha
					% vertical abaixo do ponto interno.
l=1;
while l<= npc; % for over the lines that defines the boundary
   x1 = xl1(l); y1 = yl1(l); 
   x2 = xl2(l); y2 = yl2(l);
   if(raio(l)==0) % The segment is a straight line
      if (xpi > x1 && xpi < x2) || (xpi > x2 && xpi < x1) % Check only if the x coordinate of the point is
                                        % between the x coordinate of the intial and final points of the recta
         m = (y2-y1)/(x2-x1);		% Angular coeficient of the recta
         yi = m*(xpi-x1)+y1;			% y coordinate in the intersection
         if ypi  >yi % If this is true, the point is above the line
            interv = interv + 1;	% Counter the intersection
         end;
      end;
   else   % The segment is an arc
      [xc,yc]=calcula_centro(x1,y1,x2,y2,raio(l)); % compute the center of the arc
      if(xpi<xc+abs(raio(l)) && xpi>xc-abs(raio(l))) % check only if the x coordinate of the point is between 
                 % the values xc-radius and xc+radius
         [teta_i,teta_f]=calcula_arco(x1,y1,x2,y2,xc,yc); % compute the arc between the line that defines the arc and
         % the horizontal direction (-pi<teta<pi)
         tetac(1)=acos((xpi-xc)/abs(raio(l))); % first intersection of the horizontal line that cross the 
         % point with the circunference (angle between the radius that cross the intersection point and
         % the horizontal direction)
         tetac(2)=-tetac(1);  % angle of the second intersection 
         y(1)=abs(raio(l))*sin(tetac(1))+yc; % y coordinate of the first intersection point
         y(2)=abs(raio(l))*sin(tetac(2))+yc; % y coordinate of the second intersection point
         for k=1:2 % check if the angles of the two intersection points are between the angle of the initial and the
            % final angles of the arc. If yes so the vertical line from the the point to the botton direction
            % intercept the arc (the counter should be increased).
            if(raio(l)>0) % the center is on the left of the arc (from the initial point to the final point)
               if(teta_f>teta_i) 
                  if(tetac(k)>teta_i && tetac(k)<teta_f)
                     if(y(k)<ypi)
                        interv=interv+1;
                     end;
                  end;
               else % teta_f<teta_i 
                  if(tetac(k)>teta_i || tetac(k)<teta_f)
                     if(y(k)<ypi)
                        interv=interv+1;
                     end;
                  end
               end;
            else % raio(l) < 0 the center is on the right of the arc (from the initial point to the final point)
               if(teta_i > teta_f)
                  if(tetac(k)>teta_f && tetac(k)<teta_i)
                     if(y(k)<ypi)
                        interv=interv+1;
                     end;
                  end;
               else % teta_i < teta_f
                  if(tetac(k)>teta_f || tetac(k)<teta_i)
                     if(y(k)<ypi)
                        interv=interv+1;
                     end;
                  end;
               end;
            end;
         end
      end;
   end;
   l=l+1;
end;


if rem(interv,2) ~= 0	% Resto da divisão de interv por 2
  ponto = 'interno';
else
  ponto = 'externo';
end;