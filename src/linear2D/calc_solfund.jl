function  calc_solfund(x,y,xd,yd,nx,ny,k)
# Evaluates the fundamental solution for the Laplace equation
 R=sqrt((x-xd)^2+(y-yd)^2); # Distance between the source and field points
 Rx=(x-xd); # x component of the distance
 Ry=(y-yd); # y component of the distance
 Tast=-1/(2*pi*k)*log(R); # Fundamental solution for the temperature
 qast=1/(2*pi)*(Rx*nx+Ry*ny)/R^2; # Fundamental solution for the flux
 return Tast,qast
end
# function  calc_solfund(x,y,xd,yd,nx,ny,k)
# Evaluates the fundamental solutions for the Helmholtz equation
#  r=sqrt((x-xd)^2+(y-yd)^2); # Raio (distância entre ponto fonte e ponto campo)
#  rx=(x-xd)/r; # Componente x do raio
#  ry=(y-yd)/r; # Componente y do raio
#  drdn=rx*nx+ry*ny;   #Componente do raio na direção normal
#  ZR=real(k*r);
#  Z=complex(0.,ZR);
#  F0C=SpecialFunctions.besselk(0,Z);
#  F1C=SpecialFunctions.besselk(1,Z);
#
#  qast=-(Z/r*drdn*F1C)/(2*pi); #Solução Fundamental da pressão acústica
#  Tast=F0C/(2*pi);    #Solução Fundamental do fluxo de pressão acústica
#  return Tast,qast
# end
