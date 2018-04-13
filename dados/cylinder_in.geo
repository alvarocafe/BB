// Gmsh project created on Fri Apr 13 19:13:12 2018
// Gmsh Geometry file for the inner cylinder problem
// Creating the points
s = 5.0;
Point(1) = {-5,0,0,s};
Point(2) = {0,0,0,s};
Point(3) = {5,0,0,s};
Point(4) = {-5,0,10,s};
Point(5) = {0,0,10,s};
Point(6) = {5,0,10,s};
// Creating the circles and lines
Circle(1) = {1,2,3};
Circle(2) = {4,5,6};
Line(3) = {1,4};
Line(4) = {3,6};

Circle(5) = {3,2,1};
Circle(6) = {6,5,4};

// Creating the surfaces
Line Loop(1) = {1, 4, -2, -3};
Ruled Surface(1) = {1};
Line Loop(2) = {-6, 3, 5, -4};
Ruled Surface(2) = {2};
Line Loop(3) = {2, 6};
Plane Surface(3) = {3};
Line Loop(4) = {-5, -1};
Plane Surface(4) = {4};
