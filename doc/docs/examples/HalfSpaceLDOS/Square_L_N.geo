////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

// side length 
DefineConstant [ L = 0.10 ];

// number of elements per side length
DefineConstant [ N = 4    ];

// coordinates of lower left corner
DefineConstant [ X0 = 0 ];
DefineConstant [ Y0 = 0 ];
DefineConstant [ Z0 = 0 ];

////////////////////////////////////////////////////
// end of user-specifiable parameters 
////////////////////////////////////////////////////
Point(1) = { X0,   Y0, Z0};
Point(2) = { X0+L, Y0, Z0};

Line(1)  = { 1, 2 };

Transfinite Line{1} = N+1;

Extrude { 0, L, 0  } { Line{1}; Layers{N}; }
