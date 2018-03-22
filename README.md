# BEM_base
---
This is a Julia language implementation of the boundary element method (BEM) using different Lagrangian elements.
The elements uses shape functions which are linear and quadratic.
The goal of the BEM_base project is to provide a platform that can be used alongside FreeCad and Gmsh to quickly solve bidimensional boundary element problems.
Initially, this implementation will solve the Helmholtz and Laplace equations.

To solve a problem with a geometry file 'file.msh' and boundary conditions in each face described by an array 'BCFace', to use the BEM_base to solve the problem, simply run:  
    `include("BEM_base.jl")`  
    `T,qT = BEM_base('file.msh',BCFace, k, 'Laplace') # For solving the Laplace equation`  
    `phi,qphi = BEM_base('file.msh',BCFace, k, 'Helmholtz') # For solving the Helmholtz equation`  
where k is the thermal condutivity in the Laplace equation and the wavenumber in the Helmholtz equation.
BEM_base will analyse the mesh and choose the solver accordingly, for linear and quadratic elements.

## Instalation
---
The project runs on Julia 0.5.2. To include the dependecies, simply clone the repository and include the files.  

    git clone https://github.com/alvarocafe/BEM_base  
And run the code given above.

## Contribute
---
This project is still in a very early stage and contributions are welcome.

## Support
---
If you are having issues, you can e-mail me at alvaro.cafe@hotmail.com

## License
---
This project is currently licensed under GNU GPL v.3
