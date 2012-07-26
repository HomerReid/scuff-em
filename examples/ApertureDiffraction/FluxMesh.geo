//
// GMSH mesh file for a simple planar screen parallel to the XY 
// plane and located at Z=2, i.e. a distance of 2 units from
// the aperture screen
//

//////////////////////////////////////////////////
// adjustable parameters 
//////////////////////////////////////////////////
L = 10;
Z = 2;

l = L/5; // meshing fineness

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Point(1) = {   L/2,   L/2,   Z,   l}; 
Point(2) = {   L/2,  -L/2,   Z,   l}; 
Point(3) = {  -L/2,  -L/2,   Z,   l}; 
Point(4) = {  -L/2,   L/2,   Z,   l}; 

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 1 };

Line Loop (1) = { 1, 2, 3, 4};

Plane Surface(1) = { 1 };
Physical Surface(1) = { 1 };


