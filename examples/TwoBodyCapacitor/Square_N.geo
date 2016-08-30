////////////////////////////////////////////////////
// usage example
// 
//  gmsh -2 -setnumber L 0.1 -setnumber N 4 Square_N.geo
// 
// to define a square mesh with side length 0.1 and 
// four triangles per side length 
////////////////////////////////////////////////////

// side length 
DefineConstant[ L = 1.0 ];

// number of elements per side length
DefineConstant[ N = 4   ];

// coordinates of lower left corner
X0 = 0;
Y0 = 0;

////////////////////////////////////////////////////
// end of user-specifiable parameters 
////////////////////////////////////////////////////
Point(1) = { 0, 0, 0};
Point(2) = { L, 0, 0};

Line(1)  = { 1, 2 };

Transfinite Line{1} = N+1;
Extrude { 0, L, 0  } { Line{1}; Layers{N}; }
//Physical Surface(1) = {6, 15, 19, 23, 27, 28};
