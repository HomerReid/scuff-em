//
// gmsh geometry for the unit cell of a planar screen with a
// square-lattice array of circular holes
//

//////////////////////////////////////////////////
// adjustable parameters 
//////////////////////////////////////////////////
R    = 0.14; // hole radius
L    = 0.75; // lattice spacing

l     = 0.1; // meshing fineness

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Point(0) = {     0,     0,   0,   l};
Point(1) = {     L,     0,   0,   l}; 
Line(1) = { 0, 1 };

// this extrusion defines four lines which we combine
// into a line loop
Extrude{0,L,0} { Line{1}; }
Line Loop (5) = { 1, 4, -2, -3 };

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Point(10) = { 0.5*L,     0.5*L,      0,   l};
Point(11) = { 0.5*L + R, 0.5*L,      0,   l};
Point(12) = { 0.5*L    , 0.5*L + R,  0,   l};
Point(13) = { 0.5*L - R, 0.5*L,      0,   l};
Point(14) = { 0.5*L    , 0.5*L - R,  0,   l};

Circle(15) = { 11, 10, 12 };
Circle(16) = { 12, 10, 13 };
Circle(17) = { 13, 10, 14 };
Circle(18) = { 14, 10, 11 };

Line Loop (6) = { 15, 16, 17, 18 };

//////////////////////////////////////////////////
// subtract the hole from the square to define the 
// unit cell geometry
//////////////////////////////////////////////////
Plane Surface(1) = { 5, 6 };

Translate { L,  L, 0 } { Duplicata { Surface{1}; } }
Translate { L,  0, 0 } { Duplicata { Surface{1}; } }
Translate { L, -L, 0 } { Duplicata { Surface{1}; } }
Translate { 0,  L, 0 } { Duplicata { Surface{1}; } }
Translate { 0, -L, 0 } { Duplicata { Surface{1}; } }
Translate { -L,  L, 0 } { Duplicata { Surface{1}; } }
Translate { -L,  0, 0 } { Duplicata { Surface{1}; } }
Translate { -L, -L, 0 } { Duplicata { Surface{1}; } }

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Point(1001) = {   2*L,   2*L,   1,   l};
Point(1002) = {   2*L,    -L,   1,   l};
Point(1003) = {    -L,    -L,   1,   l};
Point(1004) = {    -L,   2*L,   1,   l};
Line(1001) = { 1001, 1002 };
Line(1002) = { 1002, 1003 };
Line(1003) = { 1003, 1004 };
Line(1004) = { 1004, 1001 };
Line Loop(1005) = { 1001, 1002, 1003, 1004};

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Point(2001) = {     L,     L,   1,   l};
Point(2002) = {     L,     0,   1,   l};
Point(2003) = {     0,     0,   1,   l};
Point(2004) = {     0,     L,   1,   l};
Line(2001) = { 2001, 2002 };
Line(2002) = { 2002, 2003 };
Line(2003) = { 2003, 2004 };
Line(2004) = { 2004, 2001 };
Line Loop(2005) = { 2001, 2002, 2003, 2004};
