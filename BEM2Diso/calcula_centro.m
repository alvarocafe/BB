function [xc,yc]=calcula_centro(x1,y1,x2,y2,raio)

% Compute the center of an arc given two points and the radius

%Ponto médio Pm (xm,ym).
xm=(x1+x2)/2;
ym=(y1+y2)/2;
%Comprimento da corda (b).
b=sqrt((x2-x1)^2+(y2-y1)^2);
t1=(x2-x1)/b;
t2=(y2-y1)/b;
n1=t2;
n2=-t1;
%diferença diametro e corda (h).
h=sqrt(raio^2-(b/2)^2);
if(raio>0)
   if(n1==0)
      xc=xm;
      yc=ym-n2/abs(n2)*h;
   else
      xc=-n1/abs(n1)*sqrt(h^2*n1^2/(n1^2+n2^2))+xm;         
      yc=n2/n1*(xc-xm)+ym;
   end;
else 
   if(n1==0)
      xc=xm;
      yc=ym+n2/abs(n2)*h;
   else
      xc=n1/abs(n1)*sqrt(h^2*n1^2/(n1^2+n2^2))+xm;            
      yc=n2/n1*(xc-xm)+ym;
   end;
end;