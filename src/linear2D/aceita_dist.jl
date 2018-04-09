% Função para testar a proximidade de um ponto interno ao contorno
% Se o ponto estiver mais proximo que o desejado (d<d_min) o ponto
% é reprovado (aceita = não). Caso contrrário (aceita = sim)
%
%   Autor: Frederico Lourenço
%   Data de criação setembro de 1999 
%   Revisão 0.0

function [aceita] = aceita_dist(xpi,ypi,xl1,yl1,xl2,yl2,d_min,raio)

pc = 0;					% Contador de pontos do contorno
npc = length(xl1);	% N. de pontos do contorno
aceita = 'sim';

% Verificando a proximidade do ponto interno aos pontos do contorno		
while (strcmp(aceita,'sim') && pc < npc)
  pc = pc + 1;
  x1 = xl1(pc);
  y1 = yl1(pc);
  d = sqrt((xpi-x1)^2+(ypi-y1)^2);
  if d < d_min
    aceita = 'não';
  end;
end;

% Verificando a proximidade às linhas do contorno
l = 0;	% Contador das linhas do contorno

while (strcmp(aceita,'sim') && l < npc)
   l = l + 1;
   x1 = xl1(l); y1 = yl1(l);
   x2 = xl2(l); y2 = yl2(l);   
   if(raio(l)==0) % The segment is a straight line
     
      if (xpi > x1 && xpi < x2) || (xpi > x2 && xpi < x1)
         m = (y2-y1)/(x2-x1);		% Angular coeficient of the recta
         yi = m*(xpi-x1)+y1;			% y coordinate of the intersection of a recta normal to the current boundary
         % recta that cross the point
         dy = ypi - yi; 
         d = abs(dy*cos(atan(m))); % distance from the point to the recta
         if d < d_min 
            aceita = 'não';    
         end;
      end;
     
      if (ypi > y1 && ypi < y2) || (ypi > y2 && ypi < y1)
         if x1 == x2
            d = abs(xpi-x1);
         else 
            m = (y2-y1)/(x2-x1);			% Angular coeficient of the recta
            xi = 1/m*(ypi-y1)+x1;	% x coordinate of the intersection of a recta normal to the current boundary
    % recta that cross the point            
            dx = xpi - xi;
            d = abs(dx*sin(atan(m))); % distance from the point to the recta
         end;
         if d < d_min
            aceita = 'não';    
         end;
      end;
   else    % The segment is an arc
      [xc,yc]=calcula_centro(x1,y1,x2,y2,raio(l)); % Center of the arc
      [teta_i,teta_f]=calcula_arco(x1,y1,x2,y2,xc,yc); % angle of the lines that defines the arc with the horizontal direction
      [teta_i,teta_p]=calcula_arco(x1,y1,xpi,ypi,xc,yc);  % teta_p angle of the line that cross the center point and the
      % internal point with the horizontal direction
      if(raio(l)>0) % The center is in the left side of the arc (from the initial to the end point)
         if(teta_f>teta_i) 
            if(teta_p>teta_i && teta_p<teta_f)
               d=abs(raio(l)-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
               if d < d_min
                  aceita = 'não';    
               end;
            end;
         else % teta_f<teta_i
            if(teta_p>teta_i || teta_p<teta_f)
               d=abs(raio(l)-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
               if d < d_min
                  aceita = 'não';    
               end;
            end;
         end;
      else % raio(l) < 0 % The center is in the right side of the arc (from the initial to the end point)
         if(teta_i > teta_f)
            if(teta_p>teta_f && teta_p<teta_i)
               d=abs(abs(raio(l))-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
               if d < d_min
                  aceita = 'não';    
               end;
            end;
         else % teta_i < teta_f
            if(teta_p>teta_f || teta_p<teta_i)
               d=abs(abs(raio(l))-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
               if d < d_min
                  aceita = 'não';    
               end;
            end;
         end;
      end;         
   end;
end;


% Retorna a aprovação ou não do ponto interno