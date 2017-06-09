function [tet1,tet2] = calcula_arco(x1,y1,x2,y2,xc,yc,raio)
% Function to compute the tet1 angle between the line from point (x1,y1) 
% to (xc,yc) and the horizontal direction and the the tet2 angle between 
% the line from point (x2,y2) to (xc,yc) and the horizontal direction


dx1 = x1 - xc; dy1 = y1 - yc;
dx2 = x2 - xc; dy2 = y2 - yc;
tet1=atan2(dy1,dx1);
a=[dx1,dy1,0];
b=[dx2,dy2,0];
angle = atan2(norm(cross(a,b)),dot(a,b));
if(raio>0)
    tet2=tet1+angle;
else
    tet2=tet1-angle;
end
return
