//
// gmsh geometry specification for regular octahedron
//
// homer reid 8/2015
// 
// The default value of the side length is L=2.07124571073112,
// which is chosen to ensure that the volume of the octahedron
// matches the volume of a sphere of radius R=1.
//
// This default may be overridden by using the following
// option on the GMSH command line:
//
//  --SetNumber L 1

//////////////////////////////////////////////////
// geometric parameters
//////////////////////////////////////////////////
DefineConstant[ L = 2.07124571073112 ]; // edge length

//////////////////////////////////////////////////
// meshing length
//////////////////////////////////////////////////
l = 0.35;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
LOR2 = L/Sqrt(2.0);
Point(1) = { LOR2,      0,     0, l};
Point(2) = {    0,   LOR2,     0, l};
Point(3) = {    0,      0, +LOR2, l};
Point(4) = {    0,      0, -LOR2, l};

Line(11) = { 1, 2 };
Line(12) = { 2, 3 };
Line(13) = { 3, 1 };
Line Loop(11) = {11,12,13};
Ruled Surface(11) = {11};

Line(21) = { 1, 2 };
Line(22) = { 2, 4 };
Line(23) = { 4, 1 };
Line Loop(21) = {21,22,23};
Ruled Surface(21) = {21};

Rotate{ {0,0,1}, {0,0,0}, 1*Pi/2 } { Duplicata{ Surface{11,21};} }
Rotate{ {0,0,1}, {0,0,0}, 2*Pi/2 } { Duplicata{ Surface{11,21};} }
Rotate{ {0,0,1}, {0,0,0}, 3*Pi/2 } { Duplicata{ Surface{11,21};} }
