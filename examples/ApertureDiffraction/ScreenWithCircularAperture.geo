//
// gmsh geometry for a planar screen with a circular aperture
//
// the screen is square, with side length L, and lies in the XY plane
// centered at the origin 
//
// the circular aperture has radius R and is centered at the origin
//

//////////////////////////////////////////////////
// adjustable parameters 
//////////////////////////////////////////////////
R = 1;
L = 20;

lInner =    R/2; // meshing fineness near aperture
lOuter = 10*R/2; // meshing fineness near outer edges

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Point(1) = {   L/2,   L/2,   0,   lOuter }; 
Point(2) = {   L/2,  -L/2,   0,   lOuter }; 
Point(3) = {  -L/2,  -L/2,   0,   lOuter }; 
Point(4) = {  -L/2,   L/2,   0,   lOuter }; 

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 1 };

Line Loop (1) = { 1, 2, 3, 4};

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Point(10) = {     0,     0,   0,   lInner }; 
Point(11) = {     R,     0,   0,   lInner }; 
Point(12) = {     0,     R,   0,   lInner }; 
Point(13) = {    -R,     0,   0,   lInner }; 
Point(14) = {     0,    -R,   0,   lInner }; 

Circle(15) = { 11, 10, 12 };
Circle(16) = { 12, 10, 13 };
Circle(17) = { 13, 10, 14 };
Circle(18) = { 14, 10, 11 };

Line Loop (2) = { 15, 16, 17, 18 };

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Plane Surface(1) = { 1, 2 };
Physical Surface(1) = { 1 };


