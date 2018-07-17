# This is a function to generate plane triangular meshes for a geometry compatible with vocal tract geometries. 

function vocal_tract(Ac,file,k)
	A = Ac.*100 # From cm^2 to mm^2
	R=sqrt.(A./pi)	# Determine the radius of a circle with the corresponding area.
	filegeo = string(file,".geo")
	open(filegeo,"w") do f
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
// Definig the surface area of the glotis
Line Loop($(L)) = {-$(c+1), -$(c)};
Plane Surface($(s)) = {$(L)};
"		
		write(f,teste)
		# Update the iterators
		p +=3
		c +=4
		L +=1
		s +=1

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
	teste="
// Definig the surface area of the glotis
Line Loop($(L)) = {$(c-3), $(c-4)};
Plane Surface($(s)) = {$(L)};
"
		# Update the iterators
		L +=1
		s +=1
	write(f,teste);
	end
	println("Generating mesh...")
	run(`gmsh -2 $filegeo -v 0`)
	filemsh = string(file,".msh")
	mshinfo = const3D_tri.lermsh(filemsh,3) #Read the mesh generated
	# Build the domain points
	L = 140; # Length of the vocal tract
	n_pint = 50; # Number of domain points
	PONTOS_int = zeros(n_pint,4);
	delta = 1; # distance from both ends 
	passo = (L-2*delta)/(n_pint-1);
	for i = 1:n_pint
	    PONTOS_int[i,:] = [i 0 0 delta+(i-1)*passo];
	end
	# Set the boundary conditions for each face. Vowel /A/ model has 30 faces
	BCFace = ones(30,3);
	BCFace[:,3] = 0;
	BCFace[1,:] = [1 1 1]; # Neumann (flux = 1) to the Glotis
	BCFace[30,:] = [30 0 0]; # Dirichlet (pressure = 0) to the mouth
	u,q,uint,qint = const3D_tri.solve(mshinfo,PONTOS_int,BCFace,k)
	run(`rm $filegeo $filemsh`)
return qint
end
