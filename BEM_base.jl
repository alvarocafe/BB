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

# This is the main program of the BEM_base. It defines a function which will choose the right kernel for the element type of the mesh. 
include("src/const3D_quad/const3D_quad.jl")
include("src/const3D_tri/const3D_tri.jl")
using const3D_tri
using const3D_quad

function BEM_base(file,PONTOS_int=[],BCFace = [],k=1, equation = "wave")
    println("Importing mesh...")
    @time mshinfo = const3D_tri.lermsh(file,3) #Read the mesh generated 
    NOS_GEO,ELEM,elemint,CDC = mshinfo

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
