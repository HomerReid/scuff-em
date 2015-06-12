//////////////////////////////////////////////////
// GMSH geometry description for the unit cell of
// a silicon beam with rectangular cross section that 
// is infinitely extended in the x direction
//
// Homer Reid 20150515
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

Extrude { 0, LY, 0  } { Line{1};     Layers{NY}; }
Extrude { 0,  0, LZ } { Line{1,2};   Layers{NZ}; }
Extrude { 0, LY, 0  } { Line{6};     Layers{NY}; }
