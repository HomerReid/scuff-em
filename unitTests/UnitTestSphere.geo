//
// gmsh geometry specification for a sphere of radius 1
// used for SCUFF-EM unit tests 
//
// homer reid    8/2012
//

//************************************************************
//* input parameters      
//************************************************************
R = 0.75; // radius
l = R/2;  // meshing fineness 

//************************************************************
//* upper sphere *********************************************
//************************************************************
Point(0) = {  0 ,    0,   0,  lUpper};
Point(1) = {  0 ,    0,   R,  lUpper};
Point(2) = {  R ,    0,   0,  lMedium};
Point(3) = {  0 ,    0,  -R,  lLower};
Ellipse(1) = { 1, 0, 2, 2 };
Ellipse(2) = { 2, 0, 3, 3 };

Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{1,2}; }
Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{3,6}; }
Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{9,12}; }

Physical Surface(1) = { 5, 8, 11, 14, 17, 20 };
