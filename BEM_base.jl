# Boundary element method base (BEM_base)
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
# 

# This program follows the following scheme:
#S----P-----------------------------------M----------------------post
#^You are here-------------------------------------------------------
# This is the start of the program, S. From here, you'll need to define the
#problem P which will be solved and the matricial system M which will solve it.
#g are control points, n are the parametric
#curve generating functions, b are the boundary conditions and problem statement
#function kernel k. The 

#where S is the start of the project, G are geometry information, B are boundary conditions and problem definition (including the kernel), 
# M: 
#S----G-----------------B-----------------M----------------------post
#S----g--------n--------B-----------------M----------------------post
#S----g--------n--------b--------k--------M----------------------post
#S----P-----------------------------------c--------H-------------post
#S----P-----------------------------------c--------h--------a----post
#S----g--------n--------b--------k--------c--------h--------a----post.
# This is the main program of the BEM_base. Add here your new implementation!
# Notes:
# Returns the potential and flux on boundary points and then on domain points.
# 


include("src/const3D_quad/const3D_quad.jl")
include("src/const3D_tri/const3D_tri.jl")
using const3D_tri
using const3D_quad

function BEM_base(file,PONTOS_int=[],BCFace = [],k=1, equation = "wave")
    println("Importing mesh...")
    @time mshinfo = const3D_tri.lermsh(file,3) #Read the mesh generated 
    NOS_GEO,ELEM,elemint,CDC = mshinfo
    # Choose the right kernel for the element type of the mesh. 
    if equation == "wave"
	if size(ELEM,2) == 5
	    u,q,uint,qint = const3D_tri.solve(mshinfo,PONTOS_int,BCFace,k)
	else
	    u,q,uint,qint = const3D_quad.solve(mshinfo,PONTOS_int,BCFace,k)
	end
    elseif equation == "heat"
	if size(ELEM,2) == 5
	    u,q,uint,qint = const3D_tri_POT(mshinfo,PONTOS_int,BCFace,k)
	else
	    u,q,uint,qint = const3D_quad_POT(mshinfo,PONTOS_int,BCFace,k)
	end
    elseif typeof(equation) != String
	println("Error: the equation which will be solved must be specified.")
    end
    return u,q,uint,qint
end
