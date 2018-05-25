function calc_dsolfund(x,y,z,xd,yd,zd,n,k)

rx=x-xd;
ry=y-yd;
rz=z-zd;

r =sqrt(rx^2+ry^2+rz^2);
ZW=complex(0.,-k*r);
drdn=(rx*n[1] + ry*n[2] + rz*n[3])/r;
U=exp(ZW)/complex(4*pi*r,0.);
drdn=(rx*n[1] + ry*n[2] + rz*n[3])/r;
phiast=(ZW-complex(1.,0.))*U*complex(drdn/r,0.);
qast = (2-2*ZW-((k^2)*(r^2)))*U*complex(drdn/(r^2),0.);
return phiast,qast
end

