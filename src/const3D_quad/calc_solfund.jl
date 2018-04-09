function calc_solfund(x,y,z,xd,yd,zd,n,FR,CW)

rx=x-xd;
ry=y-yd;
rz=z-zd;

r =sqrt(rx^2+ry^2+rz^2);
ZW=complex(0.,-FR*r/CW);
U=exp(ZW)/complex(4*pi*r,0.);
drdn=(rx*n[1] + ry*n[2] + rz*n[3])/r;
Q=(ZW-complex(1.,0.))*U*complex(drdn/r,0.);
return U,Q
end

