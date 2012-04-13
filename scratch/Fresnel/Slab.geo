R = 10.0;
T = 10.0;

lf = R/4;
lc = R/1;

// 
Field[1] = MathEval;
Field[1].F = "(1 + 3*Sqrt(x*x+y*y)/10.0)*(1 - 3*z/10.0)";
Background Field = 1;

Point(1) = {0,0,0,lf};
Point(2) = {R,0,0,lc};
Point(3) = {0,R,0,lc};
Point(4) = {-R,0,0,lc};
Point(5) = {0,-R,0,lc};

Circle(1)= { 2, 1, 3 };
Circle(2)= { 3, 1, 4 };
Circle(3)= { 4, 1, 5 };
Circle(4)= { 5, 1, 2 };
Line Loop(5) = {1,2,3,4};

Plane Surface(6) = {5};
Extrude{0, 0, -T} { Surface{6}; }

Physical Surface(1) = {6, 15, 19, 23, 27, 28};
