//
// gmsh geometry for the unit cell of a planar screen with a
// square-lattice array of circular holes
//

//////////////////////////////////////////////////
// adjustable parameters 
//////////////////////////////////////////////////
R    = 0.25; // hole radius
L    = 1.00; // lattice spacing

l     = 0.05; // meshing fineness

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
Physical Surface(1) = { 1 };
