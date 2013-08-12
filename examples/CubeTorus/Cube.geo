////////////////////////////////////////////////// 
// gmsh geometry specification for cube
// homer reid
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// geometric parameters 
//////////////////////////////////////////////////
L = 0.4;   // side length

//////////////////////////////////////////////////
// this factor may be increased or decreased to   
// make the meshing more coarse or more fine in a
// uniform way over the entire object 
//////////////////////////////////////////////////
#Mesh.CharacteristicLengthFactor=1.3;

//////////////////////////////////////////////////
// these factors may be configured separately
// to make the meshing more coarse or more fine in
// particular regions of the object 
//////////////////////////////////////////////////
lCoarse =  0.05;
lFine   =  0.05;

//////////////////////////////////////////////////
// geometric description of cube /////////////////
//////////////////////////////////////////////////
Point(1) = { L/2, -L/2,  L/2, lCoarse};
Point(2) = {-L/2, -L/2,  L/2, lCoarse};
Point(3) = {-L/2, -L/2, -L/2, lFine};
Point(4) = { L/2, -L/2, -L/2, lFine};

Line(12) = {1,2};
Line(23) = {2,3};
Line(34) = {3,4};
Line(41) = {4,1};

Line Loop(1)={12, 23, 34, 41};
Plane Surface(2)={1};

out[]=Extrude{0,L,0} {Surface{2};};

Physical Surface(1)={2, 50, 54, 58, 62, 63};
