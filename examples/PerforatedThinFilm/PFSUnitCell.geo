//
// GMSH geometry file for (a simplified version of) the unit cell 
// of the structure considered by Martin-Moreno et al, PRL *86* 1114 (2001).
//
// The original structure was a silver slab of thickness 320 nm 
// with a square-lattice array of holes. This structure is the
// same lattice but on an infinitely thin PEC sheet.
//
// homer reid 10/2014
//

//////////////////////////////////////////////////
// parameters 
//////////////////////////////////////////////////
L = 0.750;    // square side length
R = 0.140;    // hole radius

l = R/2;      // meshing fineness

// square
Point(100)  = { 0, 0, 0, l};  
Point(101)  = { L, 0, 0, l};
Line(100)   = { 100, 101};
Extrude {0,L,0}  { Line{100}; }
Line Loop(100) = { 100, 103, -101, -102 };

// hole
Point(200)  = {  L/2,      L/2,   0, l };
Point(201)  = {  L/2+R,    L/2,   0, l };
Point(202)  = {    L/2,  L/2+R,   0, l };
Point(203)  = {  L/2-R,    L/2,   0, l };
Point(204)  = {    L/2,  L/2-R,   0, l };
Circle(201) = { 201, 200, 202 };
Circle(202) = { 202, 200, 203 };
Circle(203) = { 203, 200, 204 };
Circle(204) = { 204, 200, 201 };
Line Loop(200) = { 201, 202, 203, 204};

// mesh = square - hole
Plane Surface(300) = { 100, 200 };

Physical Surface(1) = {300};
