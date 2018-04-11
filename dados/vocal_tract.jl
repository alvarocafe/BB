# Julia script for generating a geo file for Gmsh
# Author: Álvaro Campos Ferreira
# Vocal tract model with adaptive meshing
# First geometry = /a/ vowel
Ac = [1.8 1.8 2.8 1.5 0.8 1.3 1.5 1.7 2.8 4.5 7.1 9.3 13.5 15.5 7.8] # From Clement (2007), the area of successive slices of the vocal tract during the phonation of the /a/ vowel by a subject, as determined by a MRI scan analysis.
A = Ac.*100 # From cm^2 to mm^2
R=sqrt.(A./pi)	# Determine the radius of a circle with the corresponding area.

# Writing the geo file

open("teste.geo","w") do f
# Initialize the iterators
p = 1	# Point iterator
c = 1	# Circle iterator
L = 1	# Line Loop iterator
s = 1	# Surface iterator
face = 10
	for i = 1:size(R,2)
		teste = "
// Defining the points
Point($(p)) = {-$(R[i]), 0, $((i-1)*face), 10.0};
Point($(p+1)) = {0, 0, $((i-1)*face),10.0};
Point($(p+2)) = {$(R[i]),0,$((i-1)*face),10.0};
// Defining the circles
Circle($(c)) = {$(p), $(p+1), $(p+2)};
Circle($(c+1)) = {$(p+2), $(p+1), $(p)};
// Defining the surface
Line Loop($(L)) = {$(c+1), $(c)};
		"
		# Update the iterators
		p +=3
		c +=2
		L +=1

		write(f,teste)
	end

end
