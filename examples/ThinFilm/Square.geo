//
// GMSH geometry file for a planar square
//
// homer reid 7/2012
//

L = 1;

l = L/5;       // meshing fineness

Point(100)  = { 0, 0, 0, l};  
Point(101)  = { L, 0, 0, l};
Line(100)   = { 100, 101};
Extrude {0,L,0}  { Line{100}; }

Line Loop(200) = { 100, 103, -101, -102 };

Plane Surface(300) = { 200 };

Physical Surface(1) = { 300 };
