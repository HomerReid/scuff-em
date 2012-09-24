//
// GMSH geometry file for the unit cell of the structure considered
// by Martin-Moreno et al,  PRL *86* 1114 (2001).
//
// This is a 750nm x 750nm square region of a slab of thickness 320 nm, 
// with a hole of diameter 280 nm etched through the center of the 
// square.
//
// homer reid 7/2012
//

//////////////////////////////////////////////////
// parameters 
//////////////////////////////////////////////////
L = 0.750;    // square side length
T = 0.320;    // slab thickness 
R = 0.140;    // hole radius

l = R/2;      // meshing fineness

//////////////////////////////////////////////////
// geometry description
//////////////////////////////////////////////////

// first create a square by extruding a line. this ensures that, 
// when we mesh the square into triangular panels, each panel edge 
// on the left and lower boundary will have a partner edge on the 
// right and upper boundary.
Point(100)  = { 0, 0, 0, l};  
Point(101)  = { L, 0, 0, l};
Line(100)   = { 100, 101};
Extrude {0,L,0}  { Line{100}; }

// the extrusion creates lines 101, 102, 103, from which
// we create a line loop as follows. (note: to my knowledge,
// the only way to discover the numbers and directions that
// GMSH assigns to the extruded lines is to open the .geo
// file in the GMSH gui and read off the line numbers from the 
// display.
Line Loop(100) = { 100, 103, -101, -102 };

// now create the hole
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

// now subtract the circle from the hole to define a surface
Plane Surface(300) = { 100, 200 };

// extrude the surface upward to create a 3D slab
Extrude{0,0,T} { Surface{300}; }

// finally, select the surfaces that we want to include in the 
// final mesh. again, to my knowledge the only way to learn
// the numbers that GMSH assigns to the extruded surfaces is 
// to open the .geo file in the GMSH gui.
Physical Surface(1) = { 300, 342, 329, 333, 337, 341 };
