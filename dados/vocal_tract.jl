# Julia script for generating a geo file for Gmsh
# Author: Álvaro Campos Ferreira
# Vocal tract model with adaptive meshing
# First geometry = /a/ vowel
Ac = [1.8 1.8 2.8 1.5 0.8 1.3 1.5 1.7 2.8 4.5 7.1 9.3 13.5 15.5 7.8] # From Clement (2007), the area of successive slices of the vocal tract during the phonation of the /a/ vowel by a subject, as determined by a MRI scan analysis.
A = Ac.*100 # From cm^2 to mm^2
R=sqrt.(A./pi)	# Determine the radius of a circle with the corresponding area.

# Writing the geo file

open("vocal_tract.geo","w") do f
# Initialize the iterators
p = 1	# Point iterator
c = 1	# Circle and line (curve) iterator
L = 1	# Line Loop iterator
s = 1	# Surface iterator
i = 1	# Cross section iterator
face = 10	# Each cross section is separated from the previous one by this distance [mm]
teste = "
// Defining the points
Point($(p)) = {-$(R[i]), 0, $((i-1)*face), 10.0};
Point($(p+1)) = {0, 0, $((i-1)*face),10.0};
Point($(p+2)) = {$(R[i]),0,$((i-1)*face),10.0};
// Defining the circles
Circle($(c)) = {$(p), $(p+1), $(p+2)};
Circle($(c+1)) = {$(p+2), $(p+1), $(p)};
"		
write(f,teste)
# Update the iterators
p +=3
c +=4

	for i = 2:size(R,2)
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


// Defining the lines between circle $(c) and $(c-1)
Line($(c+2)) = {$(p-3), $(p)};
Line($(c+3)) = {$(p-1), $(p+2)};

// Defining the surfaces of the inner side of the vocal tract
Line Loop($(L+1)) = {$(c-4), $(c+3), -$(c), -$(c+2)};
Ruled Surface($(s)) = {$(L+1)};
Line Loop($(L+2)) = {$(c-3), $(c+2), -$(c+1), -$(c+3)};
Ruled Surface($(s+1)) = {$(L+2)};
"
		# Update the iterators
		p +=3
		c +=4
		L +=3
		s +=2

		write(f,teste)
	end
end
println("Generating mesh...")
run(`gmsh -2 vocal_tract.geo`)
