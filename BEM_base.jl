# Boundary Element Method Base
# Author: Álvaro Campos Ferreira
 # Copyright (C) 2018 Álvaro Campos Ferreira
 # 
 # This program is free software; you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation; either version 3 of the License, or (at
 # your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful, but
 # WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 # General Public License for more details.
 # 
 # You should have received a copy of the GNU General Public License
 # along with this program; if not, write to the Free Software
 # Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

# This is the main program of the BEM_base.

function BEM_base(file,BCFace = [], k = 1, equation = "wave")
println("Importing mesh...")
NOS_GEO,ELEM,elemint,CDC = lermsh(file,3) #Read the mesh generated 

if equation == "wave"
	u,q,uint,qint = const3D_tri(file,BCFace,k)
elseif equation == "heat"
	u,q,uint,qint = const3D_tri_POT(file,BCFace,k)
elseif typeof(equation) != String
	println("Error: the equation which will be solved must be specified.")
end
return u,q,uint,qint
end
