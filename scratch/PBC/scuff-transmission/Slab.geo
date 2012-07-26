//
// GMSH geometry file for a 1x1 planar square 
//
// homer reid 7/2012
//

//////////////////////////////////////////////////
// parameters 
//////////////////////////////////////////////////
L = 1;
T = 1;

l = L/5;      // meshing fineness

Point(100)  = { 0, 0, 0, l};  
Point(101)  = { L, 0, 0, l};
Line(100)   = { 100, 101};
Extrude {0,L,0}  { Line{100}; }

Line Loop(200) = { 100, 103, -101, -102};

Plane Surface(300) = { 200 };

Extrude{0,0,-T} { Surface{300}; }

Physical Surface(1) = { 104, 322};

