////////////////////////////////////////////////// 
// gmsh geometry specification for torus
// homer reid
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// geometric parameters 
//////////////////////////////////////////////////
R = 1;          // outer radius
Rho = 0.35;     // inner radius

//////////////////////////////////////////////////
// this factor may be increased or decreased to   
// make the meshing more coarse or more fine in a
// uniform way over the entire object 
//////////////////////////////////////////////////
#Mesh.CharacteristicLengthFactor=0.6;

//////////////////////////////////////////////////
// these three factors may be configured separately
// to make the meshing more coarse or more fine in
// particular regions of the object 
//////////////////////////////////////////////////
lCoarse=0.3;
lMiddle=0.3;
lFine=0.3;

//////////////////////////////////////////////////
// geometric description of torus 
//////////////////////////////////////////////////
Point(1) = {0, R,        0,lMiddle};
Point(2) = {0, R+Rho,    0,lMiddle};
Point(3) = {0, R,     -Rho,lFine};
Point(4) = {0, R-Rho,    0,lMiddle};
Point(5) = {0, R,     +Rho,lCoarse};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Extrude Surface {6, {0,0,1}, {0,0,0}, 2*Pi/3};
Extrude Surface {28, {0,0,1}, {0,0,0}, 2*Pi/3};
Extrude Surface {50, {0,0,1}, {0,0,0}, 2*Pi/3};

Physical Surface(1) = { 15, 19, 23, 27, 37, 41, 45, 49, 59, 63, 67, 71};
