// Gmsh project created on Mon Oct 07 10:57:37 2024
lc = 1;
Point(1) = {0, 0, 0, lc};
Point(2) = {.1, 0,  0, lc};
Point(3) = {.1, .3, 0, lc};
Point(4) = {0,  .3, 0, lc};
Line(1) = {1, 2};
Line(2) = {3, 2};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1,-2,3,4};
Plane Surface(1) = {1};
Physical Curve(5) = {1, 2, 4};
Physical Surface("My surface") = {1};


Point(5) = {0,.4,0,lc};
Line(5)= {4,5};
Translate {-0.02, 0, 0} { Point{5}; }
Rotate {{0,0,1}, {0,0.3,0}, -Pi/4} { Point{5}; }
Translate {0, 0.05, 0} { Duplicata{ Point{3}; } }
Line(6)={5,6};
Line(7)={6,3};
Curve Loop(2)={5,6,7,3};
Plane Surface(2)={2};
surface_b[] = Translate {0.12, 0, 0} { Duplicata{ Surface{1, 2}; } };
Extrude {0,0,.5} { Surface{surface_b[1]};}
//h=0.2;
//Extrude {0,0,h}{
//  Surface{1}; Layers{{8,2},{0.5,1} };
//}
//Extrude { {0,1,0} , {-0.1,0,0.1} , -Pi/2 } {
// Surface{29}; Layers{7}; Recombine;
//}




//Extrude {0, 0, 3} {
//  Surface{2}; 
//};

//MeshSize {14,5,6,10}= lc*4 ;






