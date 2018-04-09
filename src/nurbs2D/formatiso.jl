function calcula_arcoiso(x1,y1,x2,y2,xc,yc)
# Function to compute the tet1 angle between the line from point (x1,y1) to (xc,yc) and the
# horizontal direction and the the tet2 angle between the line from point (x2,y2) to (xc,yc)
# and the horizontal direction


dx1 = x1 - xc; dy1 = y1 - yc
dx2 = x2 - xc; dy2 = y2 - yc

# Computation of tet1
      if dy1 == 0				# The point 1 and the center have the same y coordinate
        if x1 > xc
          tet1 = 0
        else  # (x1 < xc)
          tet1 = π
        end
      elseif dx1 == 0				# The point 1 and the center have the same x coordinate
        if y1 > yc
          tet1 = π/2
        else  # (y1 < yc)
          tet1 = -π/2
        end
      else  # (dx1~=0 e dy1~=0)
        tet1 = atan(dy1/dx1);
        if dx1<0 && tet1<0
          tet1 = π + tet1
        elseif dx1 < 0 && tet1>0
          tet1 = -π + tet1
        end
      end

# Computation of tet2
      if dy2 == 0				# The point 2 and the center have the same y coordinate
        if x2 > xc
          tet2 = 0
        else  # (x2 < xc)
          tet2 = π
        end
      elseif dx2 == 0				# The point 2 and the center have the same x coordinate
        if y2 > yc
          tet2 = π/2
        else  # (y2 < yc)
          tet2 = -π/2
        end
      else  # (dx2~=0 e dy2~=0)
        tet2 = atan(dy2/dx2);
        if dx2<0 && tet2<0
          tet2 = π + tet2
        elseif dx2 < 0 && tet2>0
          tet2 = -π + tet2
        end
      end
[tet1,tet2]
end


function calcula_centroiso(x1,y1,x2,y2,raio)
# Compute the center of an arc given two points and the radius

xm=(x1+x2)/2
ym=(y1+y2)/2
b=√((x2-x1)^2+(y2-y1)^2)
t1=(x2-x1)/b
t2=(y2-y1)/b
n1=t2
n2=-t1
h=√(abs(raio^2-(b/2)^2))
if(raio>0)
   if(n1==0)
      xc=xm
      yc=ym-n2/abs(n2)*h
   else
      xc=-n1/abs(n1)*√(h^2*n1^2/(n1^2+n2^2))+xm;
      yc=n2/n1*(xc-xm)+ym
   end
else
   if(n1==0)
      xc=xm
      yc=ym+n2/abs(n2)*h
   else
      xc=n1/abs(n1)*√(h^2*n1^2/(n1^2+n2^2))+xm;
      yc=n2/n1*(xc-xm)+ym
   end
end
[xc,yc]
end


function format_dad_iso(PONTOS,SEGMENTOS,MALHA)

# Programa para formata��o dos dados de entrada
# [NOS_GEO,NOS,ELEM,tipoCDC,valorCDC,normal]
#   Autor: Lucas Campos
#include("dad_0.jl")
num_elementos=sum(MALHA[:,2])
NOS_GEO=zeros(num_elementos,2)
ELEM=zeros(UInt32,num_elementos,2)
NOS=zeros(num_elementos,2)
valorCDC=zeros(num_elementos)
tipoCDC=falses(num_elementos)
normal=zeros(2,num_elementos)
cont_nos = 0;  # Counter to the physical nodes
cont_el = 0;	# Counter to the elements (the number of physical and geometric elements is the same).
num_lin = length(SEGMENTOS[:,1]);	# N�mero de linhas no contorno
p_ini = SEGMENTOS[1,2];
crv=Array{Curve}(num_lin)

#______________________________________________________________________
# Definition of the biggest dimension of the problem
max_dl = 0
for lin = 1 : num_lin
    p1 = convert(Int32,SEGMENTOS[lin,2]);
    p2 = convert(Int32,SEGMENTOS[lin,3]);
    xp1 = PONTOS[p1,2]
    yp1 = PONTOS[p1,3]
    xp2 = PONTOS[p2,2]
    yp2 = PONTOS[p2,3]
    dl = √((xp1-xp2)^2+(yp1-yp2)^2)
    if dl > max_dl
        max_dl = dl
    end
end
#_____________________________________________________________________

no_ini=1
t=1
p2=0
no1_prox=0
while(t<=num_lin)  	# While over all lines
        # Coordinates of the initial and final PONTOS of each line
        # [x1l,y1l,x2l,y2l)
        p1  = convert(Int32,SEGMENTOS[t,2])
        p2  = convert(Int32,SEGMENTOS[t,3])
        x1l = PONTOS[p1,2]
        y1l = PONTOS[p1,3]
        x2l = PONTOS[p2,2]
        y2l = PONTOS[p2,3]
        if(SEGMENTOS[t,4]==0) # The segment is a straight line
            # Increment in x and y direction
            delta_x = x2l - x1l
            delta_y = y2l - y1l
            crv[t]= nrbmak([x1l x2l;y1l y2l; 0 0], [0, 0, 1,1]);
          else #The segment is an arc
              # Compute the center of the arc and its coordinates
              r = SEGMENTOS[t,4]
              xc,yc=calcula_centroiso(x1l,y1l,x2l,y2l,r)
              # Distance between p1 and c (r1) and between p2 and c (r2)
              r1 = √((x1l-xc)^2+(y1l-yc)^2)
              r2 = √((x2l-xc)^2+(y2l-yc)^2)
              if abs(r1-r2)<.00001*max_dl
                  # Compute the angle between the lines from point c to p1 [tet1) and c to p2 (tet2]
                  tet1,tet2 = calcula_arcoiso(x1l,y1l,x2l,y2l,xc,yc)
                  if tet2 < tet1
                      tet2 = tet2 + 2*π
                  end

                  # Angle of the sector defined by the arc
                  if SEGMENTOS[t,4] > 0
                      tet = abs(tet2-tet1)
                      sig = 1.
                  else
                      tet = 2*π-abs(tet2-tet1)
                      sig = -1.
                  end
                  delta=1e-12
                  if abs(tet)-delta <= pi/2
                    narcs = 1;                # number of arc segments
                    knots = [0, 0, 0, 1, 1, 1];
                  elseif abs(tet)-delta <= pi
                    narcs = 2;
                    knots = [0, 0, 0, 0.5, 0.5, 1, 1, 1];
                  elseif abs(tet)-delta <= 3*pi/2
                    narcs = 3;
                    knots = [0, 0, 0, 1/3, 1/3, 2/3, 2/3, 1, 1, 1];
                  else
                    narcs = 4;
                    knots = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1];
                  end

                  dtet = tet/(2*narcs);     #arc segment tet angle/2

                  # determine middle control point and weight
                  wm = cos(dtet);
                  x  = r*wm;
                  y  = r*sin(dtet);
                  xm = x+y*tan(dtet);

                  # arc segment control points
                  ctrlpt = [ x wm*xm x;    # w*x - coordinate
                            -y 0     y;    # w*y - coordinate
                             1 wm    1];   # w   - coordinate

                  # build up complete arc from rotated segments
                  coefs = zeros(3,2*narcs+1);   # nurbs control points of arc
                  angle=tet1 + dtet
                  sn = sin(angle);
                  cn = cos(angle);
                  coefs[:,1:3] = [cn -sn 0 ; sn cn 0 ; 0 0 1]*ctrlpt;     # rotate to start angle
                  angle=2*dtet
                  sn = sin(angle);
                  cn = cos(angle);
                  for n = 2:narcs
                     m = 2*n+[0, 1];
                     coefs[:,m] = [cn -sn 0 ; sn cn 0 ;0 0 1]*coefs[:,m-2];
                  end

                  # vectrans arc if necessary
                  coefs = [1 0 xc; 0 1 yc; 0 0 1]*coefs;
				  linha,coluna=size(coefs)
				  coefsz=[coefs[1:2,:];zeros(coluna)';coefs[3,:]']
                  crv[t] = nrbmak(coefsz,knots)
              else
                  error("Error in the data input file: Wrong central point")
              end
          end
        t=t+1
end                                  # end of while(t<num_lin)
return crv
end


function monta_Teqiso(tipoCDC,valorCDC,x)
# Separa fluxo e temperatura

# ncdc = n�mero de linhas da matriz CDC
# T = vetor que cont�m as temperaturas nos n�s
# q = vetor que cont�m o fluxo nos n�s

ncdc = length(tipoCDC)
T=zeros(ncdc)
q=zeros(ncdc)
for i=1:ncdc # La�o sobre as condi��es de contorno
    if tipoCDC[i] == 1 # Fluxo � conhecido
        T[i] = x[i]; # A temperatura � o valor calculado
        q[i] = valorCDC[i]; # O fluxo � a condi�ao de contorno
    else # A temperatura � conhecida
        T[i] = valorCDC[i]; # A temperatura � a condi�ao de contorno
        q[i] = x[i]; # O fluxo � o valor calculado
    end
end
return T,q
end
function mostra_geo(crvs)
  p=plot(legend=:none,aspect_ratio=:equal)
  for i in crvs
    p=nrbplot(i)
  end
  p
end
