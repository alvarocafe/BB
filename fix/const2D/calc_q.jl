function calc_q(xd,yd,fc,k)
# Evaluates the influence from concentrated acoustic sources
  n_inc = size(fc,1);	# Number of acoustic sources
  q = complex(0,0);	# Allocates the influence of the concentrated sources
  for i = 1:n_inc # Loop over the concentrated sources
    x = fc[i,2];  # x coordinate of the source
    y = fc[i,3];  # y coordinate of the source
    # Evaluate the fundamental solution
    r=sqrt((x-xd)^2+(y-yd)^2); # Distance between the source and the field points
    rx=(x-xd)/r; # x component of the distance
    ry=(y-yd)/r; # y component of the distance
    ZR=real(k*r);	# wave number for the distance between the points
    Z=complex(0.,ZR);	# let the wave number be a purely imaginary number
    F0C=SpecialFunctions.besselk(0,Z); 
    Tast=F0C/(2*pi);    # Evaluates the fundamental solution for the flux
    q += fc[i,4]*Tast;	# Evaluates the influence from the concentrated source
  end

return q
end
