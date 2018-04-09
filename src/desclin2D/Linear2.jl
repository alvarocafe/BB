# Introdução ao Método dos Elementos de Contorno
# Programa de elementos de contorno aplicado a problemas de condução de
# calor. Pode haver fontes de calor concentradas
# Tipo de elementos: lineares descontínuos
# Por: Éder Lima de Albuquerque e Álvaro Campos Ferreira
# Brasília, outubro de 2017
using SpecialFunctions
using PyPlot
include("BEM_fix.jl")
# include("dad.jl"); # Arquivo de entrada de dados
include("dad_2.jl")
i=6
k=1
 # PONTOS, SEGMENTOS, MALHA, CCSeg, PONTOS_int, FR, CW,fc,finc,phi_analytical = dad_1(i,k)
 # println("Número de elementos por  meio círculo = ",i)

NOS,NOS_const,ELEM,CDC=format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg);# formata os dados
ne = size(ELEM,1);
N11,N21=calc_fforma(-1/2); # Calcula as funcoes de forma
N12,N22=calc_fforma(1/2); # Calcula as funcoes de forma
NOS_d = zeros(2*ne,3);
ELEM_d = zeros(Int,ne,3);
CDC_d = zeros(2*ne,3);
j = 0;
for i=1:ne
    no1 = ELEM[i,2];
    no2 = ELEM[i,3]
    x1 = NOS[no1,2]
    y1 = NOS[no1,3];
    x2 = NOS[no2,2];
    y2 = NOS[no2,3];
    j=j+1;
    NOS_d[j,:] = [j N11*x1+N21*x2 N11*y1+N21*y2];
    j=j+1;
    NOS_d[j,:] = [j N12*x1+N22*x2 N12*y1+N22*y2];
    ELEM_d[i,:] = [i j-1 j];
end

#println("1. Gerando pontos internos");
#PONTOS_INT=gera_p_in(NPX,NPY,PONTOS,SEGMENTOS); # gera os pontos internos
println("2. Construindo as matrizes H e G");
npg=16;
qsi,w=Gauss_Legendre(-1,1,npg);
G,H,q = monta_GeH(ELEM_d,NOS_d,CDC,k,fc,qsi,w);
# println(H)
# println(G)

A,b=aplica_CDC(G,H,CDC);
println("4. Resolvendo o sistema linear");
x=A\(b); # Calcula as variaveis desconhecidas
println("5. Ordenando fluxo e temperatura");
T,q = Monta_Teq(x,CDC);
#Ta = [0 0 0.25 0.75 1 1 0.75 0.25];
#qa = [1 1 0 0 -1 -1 0 0];
#res = H*Ta' - G*qa';
#println("6. Calculando a temperatura nos pontos internos");
#T_pint=calc_T_pint(ELEM,NOS,PONTOS_INT,T,q_vet,k,fc); # Calcula
# a temperatura nos pontos internos
#[dTdx,dTdy]=calc_dT_pint(ELEM,NOS,PONTOS_INT,T,q_vet,k,fc); # Calcula as
#  derivadas da temperatura nos
#  pontos internos
#println("7. Criando o mapa de cor");
#mapa_de_cor(T,dTdx,dTdy,NOS,ELEM,PONTOS_INT,PONTOS,SEGMENTOS,T_pint) # Faz
#um mapa de cor
