
R = 2;    // radius of loop
T = 0.2;  // wire thickness
G = 20;   // angular extent of gap, in degrees (out of 360)

l = R/5;

GR = G * Pi / 180.0;  // GR in radians

//Point(1) = {  R*Sin(GR/2), R - R*Cos(GR/2), 0.0, l };
//Point(2) = {  R*Sin(GR/2), R - R*Cos(GR/2),  10, l };
//Line(1)  = { 1, 2 };

Point(10) = {  0, 0, 0, l };
Point(11) = {  1, 1, 1, l };
Line(1)   = { 10, 11 };

//Extrude { {0,0,1}, {0, R, 0}, 2*Pi-GR } { Line{1}; }
