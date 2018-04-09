function calc_q(xd,yd,fc,k)
# Evaluates the influence from concentrated sources
valor_fonte=fc[1];
x_fonte=fc[2];
y_fonte=fc[3];
R=sqrt((x_fonte-xd)^2+(y_fonte-yd)^2); # Distance between the source and field poitns
Tast=-1/(2*pi*k)*log(R); # Fundamental solution for the temperature
q=valor_fonte*Tast;
return q
end
