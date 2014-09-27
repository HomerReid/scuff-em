// GMSH geometry for rectangular box of dimensions (Lx, Ly, Lz) with
// one corner at the origin of coordinates 
// Homer Reid 20140926

//////////////////////////////////////////////////
// user-specified parameters 
//////////////////////////////////////////////////
Lx=1.0;        // rectangle side lengths
Ly=2.0;
Lz=1.0;

Scale = 0.2; // meshing fineness lengthscale

//////////////////////////////////////////////////
// automatically derived parameters
//////////////////////////////////////////////////
NX=Ceil(Lx/Scale); // number of triangles per side
NY=Ceil(Ly/Scale);
NZ=Ceil(Lz/Scale);

Printf("NX=%g",NX);
Printf("NY=%g",NY);
Printf("NZ=%g",NZ);

Point(1) = { 0,   0,  0};
Point(2) = { Lx,  0,  0};
Line(1)  = { 1, 2 };
Transfinite Line{1} = NX+1;

// define the bottom (xy plane) surface of the box
//  (this command produces surface #5)
Extrude {0,Ly,0} { Line{1}; Layers{NY+1}; }

// extrude the bottom surface upward to define 
// the side walls and ceiling
Extrude {0,0,Lz} { Surface{5}; Layers{NZ+1}; }

// assign physical region tags to the various
// sides of the box
// (you have to open the .geo file in the GMSH gui to 
//  read off the surface numbers it assigns to 
//  the box sides)

Physical Surface(1)  = {5};   // bottom (z=0)  surface
Physical Surface(2)  = {27};  // top    (z=Lz) surface
Physical Surface(3)  = {14};  // rear   (x=0)  surface
Physical Surface(4)  = {22};  // front  (x=Lx) surface
Physical Surface(5)  = {26};  // left   (y=0)  surface
Physical Surface(6)  = {18};  // right  (y=Ly) surface
