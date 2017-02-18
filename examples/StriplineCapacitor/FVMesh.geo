////////////////////////////////////////////////////
// user-specifiable parameters 
////////////////////////////////////////////////////

// side length 
DefineConstant[ LX = 12  ];
DefineConstant[ N  = 25  ];

ZMin = -0.5;
ZMax =  1.0;
LZ=ZMax-ZMin;
NZ=N*LZ;

// number of elements per side length
DefineConstant[ N = 4 ];

NX = N*LX;
NZ = N*LZ;

Y=0.02;

////////////////////////////////////////////////////
// end of user-specifiable parameters 
////////////////////////////////////////////////////
Point(1) = { -0.5*LX, Y, ZMin};
Point(2) = {  0.5*LX, Y, ZMin};

Line(1)  = { 1, 2 };

Transfinite Line{1} = NX+1;
Extrude { 0, 0, LZ } { Line{1}; Layers{NZ}; }
