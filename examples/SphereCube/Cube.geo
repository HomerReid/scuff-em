// geometry for cube with sides of length L centered at the origin

lc = 1.1;
L = 2.0;

Point(1) = {L/2,  L/2,-L/2,lc};
Point(2) = {-L/2, L/2,-L/2,lc};
Point(3) = {-L/2,-L/2,-L/2,lc};
Point(4) = { L/2,-L/2,-L/2,lc};
Line(1)  = { 1, 2 };
Line(2)  = { 2, 3 };
Line(3)  = { 3, 4 };
Line(4)  = { 4, 1 };
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
Hello[] = Extrude { 0,0,L  } { Surface{6}; };
Physical Surface(1) = {6, 15, 19, 23, 27, 28};

// interior reference point to get normals oriented correctly
Point(100) = {0,0,0};
Physical Point(1) = {100}; 
