function dad_1(ne,FR)
#Entrada de dados para a comparação com o BEM Isogeometrico
raio = 0.5
PONTOS  = [1   -raio  0 # cilindro com centro em (0.5,0.5) de raio 0.5
          2    raio  0];
#  Segmentos que definem a geometria
#   SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final
#                                                   Raio, tipo do elemento]
#  Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
#                           inicial para o ponto final)
#                    < 0 -> O centro � � direita do segmento (do ponto
#                           inicial para o ponto final)
#                    = 0 -> O segmento � uma linha reta
#  Tipo do elemento = 1 -> Elemento quadrático contínuo
#                   = 2 -> Elemento quadrático descontínuo
#                   = 3 -> Elemento linear contínuo
SEGMENTOS = [1 1 2 -raio  #raio 0.5
             2 2 1 -raio];
# Matriz para defini��o da malha
# MALHA =[numero do segmento, numero de elementos no segmento]
MALHA = [1 ne
         2 ne];
CCSeg=[1 1 1 0
       2 1 1 0];
#FR = 1  #Frequencia analisada
CW = 343  #Velocidade de propagação da onda
#Construindo os pontos internos
n_pint = 50
PONTOS_int = zeros(n_pint,3)
dx = 1.0
dy = 0
passo = 0.05
phi_analytical = complex(zeros(size(PONTOS_int,1),1));
for i = 1:n_pint  # Para n_pint pontos internos
  PONTOS_int[i,:] = [i  dx+i*passo dy]
  bh1 = SpecialFunctions.besselh(1,2,(FR/CW)*raio);
  bh0 = SpecialFunctions.besselh(0,2,(FR/CW)*PONTOS_int[i,2]);
  phi_analytical[i,1] = (CW/FR).*(bh0./bh1);		#solucao anali­tica pela separacao de variaveis em coordenadas cilindricas da equacao de Helmholtz
end
#Incluimos as fontes concentradas e as ondas incidentes
fc = [0 3*0.5/2 -0.5/2 1];
finc = [0 -1 0 1];
return PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc,phi_analytical
end

