//
// gmsh geometry for a planar screen with a circular aperture
//
// the screen is circular, with outer radius ROut, and lies in the XY plane
// centered at the origin 
//
// the aperture is also circular, with radius RIn, and is centered at 
// the origin
//

//////////////////////////////////////////////////
// adjustable parameters 
//////////////////////////////////////////////////
RIn  = 1;
ROut = 10;

lInner =    R/3; // meshing fineness near aperture
lOuter = ROut/3; // meshing fineness near outer edges

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Point(0) = {     0,     0,   0,   lInner };
Point(1) = {  ROut,     0,   0,   lOuter }; 
Point(2) = {     0,  ROut,   0,   lOuter }; 
Point(3) = { -ROut,     0,   0,   lOuter }; 
Point(4) = {     0, -ROut,   0,   lOuter }; 

Circle(1) = { 2, 0, 1 };
Circle(2) = { 3, 0, 2 };
Circle(3) = { 4, 0, 3 };
Circle(4) = { 1, 0, 4 };

Line Loop (1) = { 1, 2, 3, 4};

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Point(11) = {   RIn,     0,   0,   lInner }; 
Point(12) = {     0,   RIn,   0,   lInner }; 
Point(13) = {  -RIn,     0,   0,   lInner }; 
Point(14) = {     0,  -RIn,   0,   lInner }; 

Circle(15) = { 11, 0, 12 };
Circle(16) = { 12, 0, 13 };
Circle(17) = { 13, 0, 14 };
Circle(18) = { 14, 0, 11 };

Line Loop (2) = { 15, 16, 17, 18 };

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Plane Surface(1) = { 1, 2 };
Physical Surface(1) = { 1 };


