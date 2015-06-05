//////////////////////////////////////////////////
// GMSH geometry description for the unit cell of
// a silicon beam with rounded sidewalls that
// is infinitely extended in the x direction
//
// Homer Reid 5/2015
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// user-specifiable parameters
//////////////////////////////////////////////////
LX  = 1.0;   // length of unit cell
LY  = 2.0;   // width of cross section
LZ  = 1.0;   // height of cross section

N  = 3;      // number of triangles per unit length

//////////////////////////////////////////////////
// end of user-specifiable parameters
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// derived parameters
//////////////////////////////////////////////////
NX = N*LX;
NY = N*LY;
NZ = N*LZ;

Point(1) = { 0,  -0.5*LY, -0.5*LZ };
Point(2) = { LX, -0.5*LY, -0.5*LZ };
Line(1)  = { 1, 2 };
Transfinite Line{1} = NX+1;

// extrusion to define bottom surface
Extrude { 0, LY, 0  } { Line{1};     Layers{NY}; }

// extrusion to define first curved sidewall
Extrude{ {-1,0,0}, {0.0, -0.5*LY, 0.0}, Pi/2} { Line{1}; Layers{3}; }
Extrude{ {-1,0,0}, {0.0, -0.5*LY, 0.0}, Pi/2} { Line{6}; Layers{3}; }

// extrusion to define top surface
Extrude { 0, LY, 0  } { Line{10};     Layers{NY}; }

// extrusion to define second curved sidewall
Extrude{ {-1,0,0}, {0.0, 0.5*LY, 0.0},  Pi/2} { Line{14}; Layers{3}; }
Extrude{ {-1,0,0}, {0.0, 0.5*LY, 0.0},  Pi/2} { Line{18}; Layers{3}; }
