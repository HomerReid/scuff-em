//////////////////////////////////////////////////
// GMSH geometry file for a linear paul surface trap
// 
// Anton Grounds 2/2/2017
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// user-definable constants //////////////////////
//////////////////////////////////////////////////

DefineConstant[ ELCNT = 28  ];  // electrode count

DefineConstant[ rfwidth  = 216  ];  // rf width
DefineConstant[ dcwidth  = 200 ];  // dc width
DefineConstant[ dclength = 1000  ];  // dc length
DefineConstant[ gndwidth = 84 ];  // gnd width
DefineConstant[ rotwidth = 20 ];  // rotation width

DefineConstant[ gndrotgap = 5  ];  // gnd - rot gap
DefineConstant[ rotrfgap  = 5 ];  // rot - rf gap
DefineConstant[ rfrotgap  = 10  ];  // rf - rot gap
DefineConstant[ rotdcgap  = 10 ];  // rot - dc gap
DefineConstant[ dcdcgap   = 20 ];  // dc - dc gap

Resgnd = gndwidth / 4; // Mesh Resolution
Resrf = rfwidth / 4;
Resrot = rotwidth / 2;
Resdc = dcwidth / 4;

MeshRes = 50;

rflength = ((dcwidth + dcdcgap) * ELCNT) - dcdcgap;

//////////////////////////////////////////////////
// Build Gnd Electrode
//////////////////////////////////////////////////
gndpoints0[] = {-rflength/2,-gndwidth/2};
gndpoints1[] = {-rflength/2,gndwidth/2};
gndpoints2[] = {rflength/2,gndwidth/2};
gndpoints3[] = {rflength/2,-gndwidth/2};

Point(1) = {gndpoints0[0],gndpoints0[1],0,Resgnd};
Point(2) = {gndpoints1[0],gndpoints1[1],0,Resgnd};
Point(3) = {gndpoints2[0],gndpoints2[1],0,Resgnd};
Point(4) = {gndpoints3[0],gndpoints3[1],0,Resgnd};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Transfinite Line {1,3} = (gndwidth/MeshRes)+1;
Transfinite Line {2,4} = (rflength/MeshRes)+1;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface("GND") = {1};
Transfinite Surface{1};
//Recombine Surface{0};


//////////////////////////////////////////////////
// Build Inner Rotation Electrodes
//////////////////////////////////////////////////

// Rotation 1

inrotpoints0[] = {gndpoints0[0],gndpoints1[1]+gndrotgap};
inrotpoints1[] = {gndpoints1[0],gndpoints1[1]+gndrotgap+rotwidth};
inrotpoints2[] = {gndpoints2[0],gndpoints1[1]+gndrotgap+rotwidth};
inrotpoints3[] = {gndpoints3[0],gndpoints1[1]+gndrotgap};

Point(5) = {inrotpoints0[0],inrotpoints0[1],0,Resrot};
Point(6) = {inrotpoints1[0],inrotpoints1[1],0,Resrot};
Point(7) = {inrotpoints2[0],inrotpoints2[1],0,Resrot};
Point(8) = {inrotpoints3[0],inrotpoints3[1],0,Resrot};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Transfinite Line {5,7} = (rotwidth/MeshRes)+1;
Transfinite Line {6,8} = (rflength/MeshRes)+1;


Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Physical Surface("Rot2") = {2};
Transfinite Surface{2};
//Recombine Surface{2};


// Rotation 2

Point(9) = {inrotpoints0[0],-inrotpoints0[1],0,Resrot};
Point(10) = {inrotpoints1[0],-inrotpoints1[1],0,Resrot};
Point(11) = {inrotpoints2[0],-inrotpoints2[1],0,Resrot};
Point(12) = {inrotpoints3[0],-inrotpoints3[1],0,Resrot};

Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,9};

Transfinite Line {9,11} = (rotwidth/MeshRes)+1;
Transfinite Line {10,12} = (rflength/MeshRes)+1;

Line Loop(3) = {9, 10, 11, 12};
Plane Surface(3) = {3};
Physical Surface("Rot3") = {3};
Transfinite Surface{3};
//Recombine Surface{3};


//////////////////////////////////////////////////
// Build RF Electrodes
//////////////////////////////////////////////////

// RF 1

rfpoints0[] = {inrotpoints0[0],inrotpoints1[1]+rotrfgap};
rfpoints1[] = {inrotpoints1[0],inrotpoints1[1]+rotrfgap+rfwidth};
rfpoints2[] = {inrotpoints2[0],inrotpoints1[1]+rotrfgap+rfwidth};
rfpoints3[] = {inrotpoints3[0],inrotpoints1[1]+rotrfgap};

Point(13) = {rfpoints0[0],rfpoints0[1],0,Resrf};
Point(14) = {rfpoints1[0],rfpoints1[1],0,Resrf};
Point(15) = {rfpoints2[0],rfpoints2[1],0,Resrf};
Point(16) = {rfpoints3[0],rfpoints3[1],0,Resrf};

Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};
Line(16) = {16,13};

Transfinite Line {13,15} = (rfwidth/MeshRes)+1;
Transfinite Line {14,16} = (rflength/MeshRes)+1;

Line Loop(4) = {13, 14, 15, 16};
Plane Surface(4) = {4};
Transfinite Surface{4};
//Recombine Surface{4};


// RF 2

Point(17) = {rfpoints0[0],-rfpoints0[1],0,Resrf};
Point(18) = {rfpoints1[0],-rfpoints1[1],0,Resrf};
Point(19) = {rfpoints2[0],-rfpoints2[1],0,Resrf};
Point(20) = {rfpoints3[0],-rfpoints3[1],0,Resrf};

Line(17) = {17,18};
Line(18) = {18,19};
Line(19) = {19,20};
Line(20) = {20,17};

Transfinite Line {17,19} = (rfwidth/MeshRes)+1;
Transfinite Line {18,20} = (rflength/MeshRes)+1;

Line Loop(5) = {17, 18, 19, 20};
Plane Surface(5) = {5};
Physical Surface("RF") = {4,5};
Transfinite Surface{5};
//Recombine Surface{5};


//////////////////////////////////////////////////
// Build Outer Rotation Electrodes
//////////////////////////////////////////////////

// Rotation 1

outrotpoints0[] = {rfpoints0[0],rfpoints1[1]+rfrotgap};
outrotpoints1[] = {rfpoints1[0],rfpoints1[1]+rfrotgap+rotwidth};
outrotpoints2[] = {rfpoints2[0],rfpoints1[1]+rfrotgap+rotwidth};
outrotpoints3[] = {rfpoints3[0],rfpoints1[1]+rfrotgap};

Point(21) = {outrotpoints0[0],outrotpoints0[1],0,Resrot};
Point(22) = {outrotpoints1[0],outrotpoints1[1],0,Resrot};
Point(23) = {outrotpoints2[0],outrotpoints2[1],0,Resrot};
Point(24) = {outrotpoints3[0],outrotpoints3[1],0,Resrot};

Line(21) = {21,22};
Line(22) = {22,23};
Line(23) = {23,24};
Line(24) = {24,21};

Transfinite Line {21,23} = (rotwidth/MeshRes)+1;
Transfinite Line {22,24} = (rflength/MeshRes)+1;

Line Loop(6) = {21, 22, 23, 24};
Plane Surface(6) = {6};
Physical Surface("Rot1") = {6};
Transfinite Surface{6};
//Recombine Surface{6};


// Rotation 2

Point(25) = {outrotpoints0[0],-outrotpoints0[1],0,Resrot};
Point(26) = {outrotpoints1[0],-outrotpoints1[1],0,Resrot};
Point(27) = {outrotpoints2[0],-outrotpoints2[1],0,Resrot};
Point(28) = {outrotpoints3[0],-outrotpoints3[1],0,Resrot};

Line(25) = {25,26};
Line(26) = {26,27};
Line(27) = {27,28};
Line(28) = {28,25};

Transfinite Line {25,27} = (rotwidth/MeshRes)+1;
Transfinite Line {26,28} = (rflength/MeshRes)+1;

Line Loop(7) = {25, 26, 27, 28};
Plane Surface(7) = {7};
Physical Surface("Rot4") = {7};
Transfinite Surface{7};
//Recombine Surface{7};


//////////////////////////////////////////////////
// Build DC Electrodes
//////////////////////////////////////////////////

tmpdcpoints0[] = {-rflength/2,outrotpoints1[1]+rotdcgap};
tmpdcpoints1[] = {-rflength/2,outrotpoints1[1]+rotdcgap+dclength};
tmpdcpoints2[] = {-rflength/2+dcwidth,outrotpoints1[1]+rotdcgap+dclength};
tmpdcpoints3[] = {-rflength/2+dcwidth,outrotpoints1[1]+rotdcgap};	

// DC 1 Upper

Point(29) = {tmpdcpoints0[0],tmpdcpoints0[1],0,Resdc};
Point(30) = {tmpdcpoints1[0],tmpdcpoints1[1],0,Resdc};
Point(31) = {tmpdcpoints2[0],tmpdcpoints2[1],0,Resdc};
Point(32) = {tmpdcpoints3[0],tmpdcpoints3[1],0,Resdc};

Line(29) = {29,30};
Line(30) = {30,31};
Line(31) = {31,32};
Line(32) = {32,29};

Transfinite Line {29,31} = (dclength/MeshRes)+1;
Transfinite Line {30,32} = (dcwidth/MeshRes)+1;

Line Loop(8) = {29, 30, 31, 32};
Plane Surface(8) = {8};
Physical Surface("UpperDC1") = {8};
Transfinite Surface{8};
//Recombine Surface{8};

// DC 1 Lower

Point(33) = {tmpdcpoints0[0],-tmpdcpoints0[1],0,Resdc};
Point(34) = {tmpdcpoints1[0],-tmpdcpoints1[1],0,Resdc};
Point(35) = {tmpdcpoints2[0],-tmpdcpoints2[1],0,Resdc};
Point(36) = {tmpdcpoints3[0],-tmpdcpoints3[1],0,Resdc};

Line(33) = {33,34};
Line(34) = {34,35};
Line(35) = {35,36};
Line(36) = {36,33};

Transfinite Line {33,35} = (dclength/MeshRes)+1;
Transfinite Line {34,36} = (dcwidth/MeshRes)+1;

Line Loop(9) = {33, 34, 35, 36};
Plane Surface(9) = {9};
Physical Surface("LowerDC1") = {9};
Transfinite Surface{9};
//Recombine Surface{9};


// Duplicate and translate
Geometry.CopyMeshingMethod = 1;
For i In {1:ELCNT-1}
	electrode = Translate {(dcwidth + dcdcgap)*i, 0, 0} { Duplicata{ Surface{8}; } };
	Physical Surface (StrCat("UpperDC", Sprintf("%g",(i+1)))) = {electrode};
	
	electrode = Translate {(dcwidth + dcdcgap)*i, 0, 0} { Duplicata{ Surface{9}; } };
	Physical Surface (StrCat("LowerDC", Sprintf("%g",(i+1)))) = {electrode};

EndFor
