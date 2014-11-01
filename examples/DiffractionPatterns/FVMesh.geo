//
// GMSH mesh for the field-visualization surface on which we
// plot transmitted field intensity.
//
L = 0.750;    // square side length

Z = 1.0;      // distance of visualization surface from screen

l = L/10;      // meshing fineness

// square
Point(100)  = { 0, 0, Z, l};  
Point(101)  = { L, 0, Z, l};
Line(100)   = { 100, 101};
Extrude {0,L,0}  { Line{100}; }
Line Loop(100) = { 100, 103, -101, -102 };

Plane Surface(1) = { 100 };
Physical Surface(1) = { 1 };
