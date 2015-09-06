//
// gmsh geometry specification for regular tetrahedron 
//
// homer reid 8/2015
// 
// The default value of the side length is L=3.2878976200992,
// which is chosen to ensure that the volume of the tetrahedron
// matches the volume of a sphere of radius R=1.
//
// This default may be overridden by using the following
// option on the GMSH command line:
//
//  --SetNumber L 1

//////////////////////////////////////////////////
// geometric parameters
//////////////////////////////////////////////////
DefineConstant[ L = 3.2878976200992 ]; // edge length
DefineConstant[ Mesh3D        = 0 ];

//////////////////////////////////////////////////
// meshing lengthscale
//////////////////////////////////////////////////
l = 0.35;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
X = L/2.0;
Y = L/(2.0*Sqrt(3.0));
Z = L*Sqrt(2.0/3.0)/4.0;

Point(1) = {   -X,       -Y,      -Z, l};
Point(2) = {    X,       -Y,      -Z, l};
Point(3) = {    0,      2*Y,      -Z, l};
Point(4) = {    0,        0,     3*Z, l};

Line(12) = { 1, 2 };
Line(13) = { 1, 3 };
Line(14) = { 1, 4 };
Line(23) = { 2, 3 };
Line(24) = { 2, 4 };
Line(34) = { 3, 4 };

Line Loop(11) = {12,23,-13};
Line Loop(21) = {12,24,-14};
Line Loop(31) = {13,34,-14};
Line Loop(41) = {23,34,-24};
Ruled Surface(11) = {11};
Ruled Surface(21) = {21};
Ruled Surface(31) = {31};
Ruled Surface(41) = {41};

If( Mesh3D==1 )
  Surface Loop(51) = { 11, 21, 31, 41 };
  Volume(1000) = { 51 };
EndIf
