//
// gmsh geometry specification for hemispherical shell
// used to visualize radiation pattern for dipole and Yagi
// antennas
// 
// homer reid 
//

//************************************************************
//* input parameters      
//************************************************************
R = 500.0;
l = R/5;   // meshing fineness

//************************************************************
//************************************************************
//************************************************************
Point(1) = {  0,  0, 0, l };
Point(2) = {  0, -R, 0, l };
Point(3) = {  R,  0, 0, l };
Point(4) = {  0,  R, 0, l };
Point(5) = { -R,  0, 0, l };

Circle(1) = { 2, 1, 3 };
Circle(2) = { 3, 1, 4 };

Circle(3) = { 4, 1, 5 };
Circle(4) = { 5, 1, 2 };

Extrude { {0,1,0}, {0,0,0}, Pi/2} { Line{1, 2}; }
Extrude { {0,1,0}, {0,0,0}, -Pi/2} { Line{1, 2}; }

Extrude { {0,1,0}, {0,0,0}, Pi/2} { Line{3,4}; }
Extrude { {0,1,0}, {0,0,0}, -Pi/2} { Line{3,4}; }

Line Loop(1) = {1,2,3,4};
Plane Surface(100) = { 1 };

//Physical Surface(1) = { 7, 10, 28, 25 };  // lower hemisphere
Physical Surface(2) = { 13, 16, 19, 22 };   // upper hemisphere
//Physical Surface(3) = { 100 };            // equatorial plane
