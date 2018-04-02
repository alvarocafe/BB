# Julia script for generating a geo file for Gmsh
# Author: Álvaro Campos Ferreira
# Vocal tract model with adaptive meshing
# First geometry = /a/ vowel
Ac = [1.8 1.8 2.8 1.5 0.8 1.3 1.5 1.7 2.8 4.5 7.1 9.3 13.5 15.5 7.8] # From Clement (2007), the area of successive slices of the vocal tract during the phonation of the /a/ vowel by a subject, as determined by a MRI scan analysis.
A = a.*100 # From cm^2 to mm^2
R=sqrt.(A./pi)	# Determine the radius of a circle with the corresponding area.

# Writing the geo file
c = 1
for i = 1:size(c,2)
# Defining the points
	println("Point($(i)) = {-$(R), 0, 0}")
	println("Point($(i+1)) = {0, 0, 0}")
	println("Point($(i+2)) = {0, 0, $(R)}")
# Defining the circles
	println("Circle($(i)) = {$(i), $(i+1), $(i+2)}")
	println("Circle($(i+1)) = {$(i+2), $(i+1), $(i)}")
# Defining the surface
	println("Line Loop($(c)) = {$(i), $(i+1), $(i+2)}")

end
