#Conjunto de funções para o Método dos Elementos de Contorno direto
using SpecialFunctions

include("format_dad.jl")
include("calcula_centro.jl")
include("calcula_arco.jl")
include("cal_HeG.jl")
include("cal_Aeb.jl")
include("calc_F.jl")
include("Gauss_Legendre.jl")
include("calc_fforma.jl")
include("aplica_CDC.jl")
include("calc_phi_pint.jl")
include("telles.jl")    # Para o tratamento da singularidade, transformada de Telles
include("nurbs.jl")  # Implementação das NURBS
include("formatiso.jl") # Suporte para gerar curvas NURBS
include("calcula_iso.jl")   # Integração de elementos isogeométricos (NURBS)
include("ACA.jl") # Ferramentas para o ACA
