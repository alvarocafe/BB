function calc_gh_nsing(x1,y1,x2,y2,xd,yd,k,qsi,w)
# Evaluates the non singular terms of matrices G and H.
npg=size(qsi,1); # Number of integration (Gaussian quadrature) points

#h=[0. 0.]; # Initializes the matrix H terms
#g=[0. 0.]; # Initializes the matrix G terms
h = complex(zeros(2)) # Initializes the matrix H terms
g = complex(zeros(2)) # Initializes the matrix G terms

d=sqrt((x2-x1)^2+(y2-y1)^2); # Distance between the two points (half the length of the element)
L = 2*d;	# Length of the element
dgamadqsi=L/2; # Jacobian

sx=(x2-x1)/d; # x component of the tangent vector
sy=(y2-y1)/d; # y component of the tangent vector
nx=sy; # x component x of the normal vector
ny=-sx; # y component of the normal vector

for kk=1:npg	# Loop over the integration points
    N1,N2=calc_fforma_d(qsi[kk]); # Evaluates the shape functions

    x=N1*x1+N2*x2; # Evaluates the x coordinate of the integration point
    y=N1*y1+N2*y2; # Evaluates the y coordinate of the integration point

    phiast,qast=calc_solfund(x,y,xd,yd,nx,ny,k); # Evaluate the fundamental solutions
    #phiast = 1;
    h=h+[N1; N2].*qast.*dgamadqsi.*w[kk]; # Integra��o da matriz h
    g=g+[N1; N2].*phiast.*dgamadqsi.*w[kk]; # Integra��o da matriz g
end
return g,h
end
