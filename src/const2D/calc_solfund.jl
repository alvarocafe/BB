function calc_solfundpot(x,y,xd,yd,nx,ny,k)
# Evaluates the fundamental solutions of the Laplace equation.

r=sqrt((x-xd)^2+(y-yd)^2); # Distance between the source and field points
rx=(x-xd); # x component of the distance
ry=(y-yd); # y component of the distance
Tast=-1/(2*pi*k)*log(r); # Fundamental solution for the temperature
qast=1/(2*pi)*(rx*nx+ry*ny)/r^2; # Fundamental solution for the flux
return Tast, qast
end



function  calc_solfund(x,y,xd,yd,nx,ny,k)
# Evaluates the fundamental solutions of the Helmholtz equation.
r=sqrt((x-xd)^2+(y-yd)^2); # Distance between the source and field points
rx=(x-xd); # x component of the distance
ry=(y-yd); # y component of the distance
drdn=rx*nx+ry*ny;   # Distance in the normal direction
ZR=real(k*r);
Z=complex(0.,ZR);
F0C=besselk(0,Z);
F1C=besselk(1,Z);

  qast=-(Z/r*drdn*F1C)/(2*pi); 	# Fundamental solution for the velocity potential
  Tast=F0C/(2*pi);    		# Fundamental solution for the flux
  return Tast,qast
end


