//
//
//

L = 1;
G = 0.2;
H = 0.1;

l = 0.1;

Point(1) = { G/2, 0, 0, l };
Point(2) = { G/2, H, 0, l };
Point(3) = { L/2, H, 0, l };
Point(4) = { L/2, L+H, 0, l };
Point(5) = { -L/2, L+H, 0, l };
Point(6) = { -L/2, H, 0, l };
Point(7) = { -G/2, H, 0, l };
Point(8) = { -G/2, 0, 0, l };

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,1};

Line Loop(1)={1,2,3,4,5,6,7,8};
Plane Surface(1) = { 1 };
