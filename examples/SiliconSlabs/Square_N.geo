////////////////////////////////////////////////////
// user-specifiable parameters 
////////////////////////////////////////////////////

// side length 
L = 2.00;

// number of elements per unit length (higher=greater resolution)
N = 3;

////////////////////////////////////////////////////
// end of user-specifiable parameters 
////////////////////////////////////////////////////

NL = Ceil(N*L);

Point(1) = { 0, 0, 0, 1};
Point(2) = { L, 0, 0, 1};

Line(1)  = { 1, 2 };

Transfinite Line{1} = NL+1;
Extrude { 0, L, 0  } { Line{1}; Layers{NL}; }
