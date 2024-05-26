// geometry-description file created 2024-05-26 14:42:08 by system76-pc
// command used:
//   domain.py -hmesh 160.0 -o coarse.geo

lc = 160.000000;
lc_corner = 40.000000;
Point(1) = {3000.000000,0.000000,0,lc};
Point(2) = {3000.000000,400.000000,0,lc};
Point(3) = {0.000000,400.000000,0,lc};
Point(4) = {0.000000,0.000000,0,lc};
Point(5) = {1500.000000,0.000000,0,lc};
Point(6) = {1500.000000,100.000000,0,lc_corner};
Point(7) = {1750.000000,100.000000,0,lc};
Point(8) = {2000.000000,100.000000,0,lc_corner};
Point(9) = {2000.000000,0.000000,0,lc};
Line(11) = {1,2};
Line(12) = {2,3};
Line(13) = {3,4};
Line(14) = {4,5};
Line(15) = {5,6};
Line(16) = {6,7};
Line(17) = {7,8};
Line(18) = {8,9};
Line(19) = {9,1};
Line Loop(21) = {11,12,13,14,15,16,17,18,19};
Plane Surface(31) = {21};
Physical Line(41) = {11};
Physical Line(42) = {12};
Physical Line(43) = {13};
Physical Line(44) = {14,15,16,17,18,19};
Physical Surface(51) = {31};
