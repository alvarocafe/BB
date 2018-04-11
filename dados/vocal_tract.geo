// Gmsh project created on Mon Apr  2 09:46:31 2018
// Author: Álvaro Campos Ferreira
// Vocal tract model with adaptive meshing
Point(1) = {0, 0, 0, 1.0};	// Create point 1 at (x,y,z) = (0,0,0)
Delete {	// Deletes point 1
  Point{1};
}
// Defining the geometry
// First circle
Point(1) = {-10, 0, 0, 10.0};
Point(2) = {0, 0, 0, 10.0};
Point(3) = {10, 0, 0, 10.0};
Circle(1) = {1, 2, 3};	// Bottom half of the circle
Circle(2) = {3, 2, 1};	// Top half of the circle
Line Loop(1) = {-2, -1};	// Select the boundary for area creation
Plane Surface(1) = {1};	// Creates the surface of the circle

// Second circle
Point(4) = {-10, 0, 100, 10.0};
Point(5) = {0, 0, 100, 10.0};
Point(6) = {10, 0, 100, 10.0};
Circle(3) = {4, 5, 6};	// Bottom half of the circle
Circle(4) = {6, 5, 4};	// Top half of the circle
Line Loop(2) = {4, 3};	// Select the boundary for area creation
Plane Surface(2) = {2};	// Creates the surface of the circle

// Create surface area of the cylinder

Line(7) = {1, 4};
Line(8) = {3, 6};

Line Loop(3) = {1, 8, -3, -7};
Ruled Surface(3) = {3};
Line Loop(4) = {2, 7, -4, -8};
Ruled Surface(4) = {4};
