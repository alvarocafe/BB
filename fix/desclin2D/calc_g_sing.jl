function calc_g_sing(x1,y1,x2,y2,k)
# Evaluates the singular term for matrix g (this subroutine is for the continuous linear element)
L = sqrt((x1-x2)^2+(y1-y2)^2);
g = L/(4*pi*k)*(3/2-log(L));
return g
