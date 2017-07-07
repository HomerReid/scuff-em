////////////////////////////////////////////////////
// GMSH geometry file for a planar rectangular screen lying in XY plane
//
// usage: gmsh -2 Square_N
//
////////////////////////////////////////////////////

// side length 
DefineConstant[ LX = 1.0 ];
DefineConstant[ LY = 1.0 ];

// center 
DefineConstant[ X0 = 0.0 ];
DefineConstant[ Y0 = 0.0 ];

// number of elements per side length
DefineConstant[ N = 4 ];

////////////////////////////////////////////////////
// end of user-specifiable parameters
////////////////////////////////////////////////////
NX = Ceil[ N*LX ];
NY = Ceil[ N*LY ];

Point(1) = { X0-0.5*LX, Y0-0.5*LY, 0};
Point(2) = { X0+0.5*LX, Y0-0.5*LY, 0};

Line(1)  = { 1, 2 };
Transfinite Line{1} = NX+1;

Extrude { 0, LY, 0  } { Line{1}; Layers{NY}; }
