# Julia script for generating a geo file for Gmsh
# Author: Álvaro Campos Ferreira
# Vocal tract model with adaptive meshing
# First geometry = /a/ vowel
include("vocal_tractfun.jl")
#Ac = [1.4 3.4 6.6 8.7 9.8 9.2 7.6 5.8 4.4 2.3 0.7 1.5 3.9 8.6 6.9 1.5 1.3]
#Ac = [1.8 1.8 2.8 1.5 0.8 1.3 1.5 1.7 2.8 4.5 7.1 9.3 13.5 15.5 7.8] # From Clement (2007), the area of successive slices of the vocal tract during the phonation of the /a/ vowel by a subject, as determined by a MRI scan analysis.
#A = Ac.*100 # From cm^2 to mm^2
#R=sqrt.(A./pi)	# Determine the radius of a circle with the corresponding area.
include("../../src/const3D_tri/const3D_tri_vc.jl")
n = 5
file = Array{Any}(n)
for i =1:n
	file[i] = "vocal_tract$(i)"
	size::Int8 = 14
	Ac = 15.*rand(size)
	vocal_tract(Ac,file[i],const3D_tri)
end




