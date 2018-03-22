## BEM_base
# Programa de Elementos de Contorno Constantes Direto 2D em Julia
# Resolve problemas descritos pela equação de Laplace e de Helmholtz
# Autor: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com

# Incluir as dependências
using SpecialFunctions  # Modulo para definir funções especiais (Bessel e Hankel)
#using PyPlot            # Modulo para produzir plots utilizando o Matplotlib
#using KrylovMethods
include("fix/BEM_fix.jl")   # Contém as funções para resolver problemas de elementos de contorno
include("dados/dad.jl") # Contem diversas entrada de dados

## Resolvendo o problema do cilindro vibrante
# O problema do cilindro vibrante consiste em um cilindro de paredes rigidas
# imerso em um fluido compressível com velocidade de propagação de onda CW.
# O cilindro é longo o suficiente para que seja possível modelar apenas uma seção
# bidimensional como um círculo em um plano. Como não há escoamento, o comportamento
# do fluido pode ser descrito pela equação da onda bidimensional. Resolvendo o problema
# no domínio da frequência e considerando o movimento harmônico com frequência FR,
# a equação da onda se reduz à equação de Helmholtz. O problema possui uma solução analítica.

# Entrada de dados
# A função dad_1 recebe o número de elementos para cada metade do cilindro e frequência do problema.
i = 100 # Numero de elementos por metade do cilindro
FR = 20 # Frequência do problema
CW = 343 # Velocidade de propagação da onda
k = FR/CW # Numero de onda, constante adimensional que relaciona a frequencia com a velocidade de propagação.
PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc,phi_analytical = dad_1(i,FR)
#Geração dos pontos e pesos de Gauss
npg=16; # Numero de pontos de integração
qsi,w = Gauss_Legendre(-1,1,npg) # Geração dos pontos e pesos de Gauss.

## Elementos Constantes
NOS_GEO,NOS,ELEM,CDC = format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg) # Constroi as matrizes de discretização
nnos = size(NOS,1)  # Numero de nos fisicos, corresponde ao numero de elementos no caso de elementos constantes
b1 = 1:nnos # Vetor com os nos e elementos de integração.
println("Construindo A e B usando o BEM tradicional.")
@time A,B = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,CDC,CW,FR,qsi,w])  # Constroi as matrizes A e B utilizando
# o método de colocação e aplica as condições de contorno.
b = B*CDC[:,3]  # Constroi o vetor b para o sistema linear
x = A\b # Resolve o sistema linear.
phi,qphi = monta_Teq(CDC,x) # Aplica as condições de contorno para obter as variáveis no contorno.
println("Obtendo os valores nos pontos internos.")
@time phi_pint = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,CW,FR,fc,finc,qsi,w) # Obtem os valores nos pontos internos
println("Calculando o erro.")
@time erro = abs((sum((phi_pint - phi_analytical).^2)))   # Calcula a norma em comparação com a solução analítica.
println("")
## Elementos lineares descontínuos
