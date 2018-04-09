function calc_F(NOS,CW,FR,fc,finc)
nnos=size(NOS,1);
q = complex(zeros(nnos,1));
inc = complex(zeros(nnos,1));
for i = 1:nnos
  xd=NOS[i,2]; # Coordenada x do ponto fonte
  yd=NOS[i,3]; # Coordenada y do ponto fonte
  q[i]=calc_q(xd,yd,fc,FR,CW);
  inc[i] = calc_inc(xd,yd,finc,FR,CW);
end
return q, inc
end

function calc_q(xd,yd,fc,FR,CW)
  n_inc = size(fc,1);
  q = complex(0,0);
  for i = 1:n_inc
    x = fc[i,2];  #Coordenadas x e y da fonte i
    y = fc[i,3];
    #Calcula as soluções fundamentais
    r=sqrt((x-xd)^2+(y-yd)^2); # Raio (distância entre ponto fonte e
    # ponto campo)
    rx=(x-xd)/r; # Componente x do raio
    ry=(y-yd)/r; # Componente y do raio
#    drdn=rx*nx+ry*ny;   #Componente do raio na direção normal
    ZR=real(FR*r/CW);
    Z=complex(0.,ZR);
    #     [F0C,F1C]=Bessel(Z);
    F0C=SpecialFunctions.besselk(0,Z);
#    F1C=besselk(1,Z);
#    qast=-(Z/r*drdn*F1C)/(2*pi); #Solução Fundamental da pressão acústica
    Tast=F0C/(2*pi);    #Solução Fundamental do fluxo de pressão acústica
    q += fc[i,4]*Tast;
  end

return q
end

function calc_inc(xd,yd,finc,FR,CW)
n_inc = size(finc,1);
inc = complex(0,0);
for i = 1:n_inc
  inc += finc[i,4]*exp(complex(0,FR/CW*(finc[i,2]*xd + finc[i,3]*yd)));
end
return inc
end
