/* lcar1 = .0275; */
lcar1 = 1.0;

length = 10.1;
height = 9.9;
depth = 10.0;

Point(newp) = {length/2,height/2,depth,lcar1}; /* Point      1 */
Point(newp) = {length/2,height/2,0,lcar1}; /* Point      2 */
Point(newp) = {-length/2,height/2,depth,lcar1}; /* Point      3 */
Point(newp) = {-length/2,-height/2,depth,lcar1}; /* Point      4 */
Point(newp) = {length/2,-height/2,depth,lcar1}; /* Point      5 */
Point(newp) = {length/2,-height/2,0,lcar1}; /* Point      6 */
Point(newp) = {-length/2,height/2,0,lcar1}; /* Point      7 */
Point(newp) = {-length/2,-height/2,0,lcar1}; /* Point      8 */
Line(1) = {3,1};
Line(2) = {3,7};
Line(3) = {7,2};
Line(4) = {2,1};
Line(5) = {1,5};
Line(6) = {5,4};
Line(7) = {4,8};
Line(8) = {8,6};
Line(9) = {6,5};
Line(10) = {6,2};
Line(11) = {3,4};
Line(12) = {8,7};
Line Loop(13) = {-6,-5,-1,11};
Plane Surface(1) = {13};
Line Loop(15) = {4,5,-9,10};
Plane Surface(2) = {15};
Line Loop(17) = {3,+12,-8,-10};
Plane Surface(3) = {17};
Line Loop(19) = {-7,-12,+2,-11};
Plane Surface(4) = {19};
Line Loop(21) = {-4,-3,-2,1};
Plane Surface(5) = {21};
Line Loop(23) = {8,9,6,7};
Plane Surface(6) = {23};

/* Surface Loop(25) = {14,24,-18,22,16,-20}; */
//+
Physical Surface("1") = {1};
//+
Physical Surface("2") = {2};
//+
Physical Surface("3") = {3};
//+
Physical Surface("4") = {4};
//+
Physical Surface("5") = {5};
//+
Physical Surface("6") = {6};
