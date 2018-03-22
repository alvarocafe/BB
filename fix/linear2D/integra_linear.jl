#Integração dos elementos lineares
#Input dos dados
# p1 = [0 0]
# p2 = [1 1]
# NOS_GEO,NOS,ELEM,CDC = dad_const(p1,p2)
PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc,phi_analytical = dad_1(4)
NOS_GEO,NOS,ELEM,CDC = format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg)
CW = 1
FR = 1
CDC = [1.0  0.0  0.0 1 1
 2.0  1.0  0.0 0 1
 3.0  1.0  0.0 0 1
 4.0  1.0  0.0 0 1
 5.0  1.0  0.0 0 1
 6.0  1.0  0.0 0 1
 7.0  1.0  0.0 0 1
 8.0  1.0  0.0 0 1]
G,H = monta_GeH_linear(ELEM,NOS_GEO,1)
A,b,T_PR = aplica_CDC(G,H,CDC,ELEM)
