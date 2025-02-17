
Fac = 1;

Lc_b=10;
Lc_b = Lc_b* Fac;

n_layer_z = 150;

x = 1000;
y = 1000;
z = 1000;
x = x* Fac;
y = y* Fac;
z = z* Fac;

stp = 0;

PML_layers = 2;

PML_w = PML_layers*x/n_layer_z;
PML_w = PML_w * Fac;

Point(01) = {-x/2 + stp,  -y/2 + stp, -z/2 -PML_w + stp,Lc_b};
Point(02) = { x/2 + stp,  -y/2 + stp, -z/2 -PML_w + stp,Lc_b};
Point(03) = { x/2 + stp,   y/2 + stp, -z/2 -PML_w + stp,Lc_b};
Point(04) = {-x/2 + stp,   y/2 + stp, -z/2 -PML_w + stp,Lc_b};

Point(05) = {-x/2 - PML_w + stp,  -y/2 - PML_w + stp, -z/2 -PML_w + stp ,Lc_b};
Point(06) = {-x/2 +     0 + stp,  -y/2 - PML_w + stp, -z/2 -PML_w + stp ,Lc_b};
Point(07) = { x/2 +     0 + stp,  -y/2 - PML_w + stp, -z/2 -PML_w + stp ,Lc_b};
Point(08) = { x/2 + PML_w + stp,  -y/2 - PML_w + stp, -z/2 -PML_w + stp ,Lc_b};

Point(09) = { x/2 + PML_w + stp ,  -y/2 -     0 + stp , -z/2 -PML_w + stp ,Lc_b};
Point(10) = { x/2 + PML_w + stp ,   y/2 -     0 + stp , -z/2 -PML_w + stp ,Lc_b};
Point(11) = { x/2 + PML_w + stp ,   y/2 + PML_w + stp , -z/2 -PML_w + stp ,Lc_b};

Point(12) = { x/2 +     0 + stp ,   y/2 + PML_w + stp , -z/2 -PML_w + stp ,Lc_b};
Point(13) = { x/2 + PML_w + stp ,   y/2 + PML_w + stp , -z/2 -PML_w + stp ,Lc_b};
Point(14) = {-x/2 +     0 + stp ,   y/2 + PML_w + stp , -z/2 -PML_w + stp ,Lc_b};
Point(15) = {-x/2 - PML_w + stp ,   y/2 + PML_w + stp , -z/2 -PML_w + stp ,Lc_b};

Point(16) = {-x/2 - PML_w + stp ,   y/2 +      0 + stp , -z/2 -PML_w + stp ,Lc_b};
Point(17) = {-x/2 - PML_w + stp ,  -y/2 +      0 + stp , -z/2 -PML_w + stp ,Lc_b};
Point(18) = {-x/2 - PML_w + stp ,  -y/2 -  PML_w + stp , -z/2 -PML_w + stp ,Lc_b};

Line(01) = {1,2};
Line(02) = {2,3};
Line(03) = {3,4};
Line(04) = {4,1};

Line(05) = {5,6};
Line(06) = {6,1};
Line(07) = {1,17};
Line(08) = {17,5};

Line(09) = {6,7};
Line(10) = {7,2};
//Line(11) = {2,1};
//Line(12) = {1,6};

Line(13) = {7,8};
Line(14) = {8,9};
Line(15) = {9,2};
//Line(16) = {2,7};

//Line(17) = {2,9};
Line(18) = {9,10};
Line(19) = {10,3};
//Line(20) = {3,2};

//Line(21) = {3,10};
Line(22) = {10,11};
Line(23) = {11,12};
Line(24) = {12,3};

//Line(25) = {4,3};
//Line(26) = {3,12};
Line(27) = {12,14};
Line(28) = {14,4};

Line(29) = {16,4};
//Line(30) = {4,14};
Line(31) = {14,15};
Line(32) = {15,16};

//Line(33) = {17,1};
//Line(34) = {1,4};
//Line(35) = {4,16};
Line(36) = {16,17};


Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Line Loop(2) = {5,6,7,8};
Plane Surface(2) = {2};

Line Loop(3) = {9,10,-1,-6};
Plane Surface(3) = {3};

Line Loop(4) = {13,14,15,-10};
Plane Surface(4) = {4};

Line Loop(5) = {-15,18,19,-2};
Plane Surface(5) = {5};

Line Loop(6) = {-19,22,23,24};
Plane Surface(6) = {6};

Line Loop(7) = {-3,-24,27,28};
Plane Surface(7) = {7};

Line Loop(8) = {29,-28,31,32};
Plane Surface(8) = {8};

Line Loop(9) = {-7,-4,-29,36};
Plane Surface(9) = {9};

Transfinite Line{10,15,19,24,28,29,7,6} = PML_layers+1;
Transfinite Line{13,14,22,23,5,8,31,32} = PML_layers+1;

Transfinite Line{27,9,1,3} = n_layer_z +1;
Transfinite Line{36,4,18,2} = n_layer_z +1;

Transfinite Surface{1}={1,2,3,4};
Transfinite Surface{2}={5,6,1,17};
Transfinite Surface{3}={6,7,2,1};
Transfinite Surface{4}={7,8,9,2};
Transfinite Surface{5}={2,9,10,3};
Transfinite Surface{6}={3,10,11,12};
Transfinite Surface{7}={4,3,12,14};
Transfinite Surface{8}={16,4,14,15};
Transfinite Surface{9}={17,1,4,16};

Recombine Surface{1};
Recombine Surface{2};
Recombine Surface{3};
Recombine Surface{4};
Recombine Surface{5};
Recombine Surface{6};
Recombine Surface{7};
Recombine Surface{8};
Recombine Surface{9};

Extrude {0,0, PML_w} {Surface{1}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{2}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{3}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{4}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{5}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{6}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{7}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{8}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{9}; Layers{PML_layers}; Recombine; }

var1 = 58;
var2 = 80;
Extrude {0,0, z} {Surface{var1}; Layers{n_layer_z}; Recombine; }
Extrude {0,0, z} {Surface{var2}; Layers{n_layer_z}; Recombine; }
Extrude {0,0, z} {Surface{var2+1*22}; Layers{n_layer_z}; Recombine; }
Extrude {0,0, z} {Surface{var2+2*22}; Layers{n_layer_z}; Recombine; }
Extrude {0,0, z} {Surface{var2+3*22}; Layers{n_layer_z}; Recombine; }
Extrude {0,0, z} {Surface{var2+4*22}; Layers{n_layer_z}; Recombine; }
Extrude {0,0, z} {Surface{var2+5*22}; Layers{n_layer_z}; Recombine; }
Extrude {0,0, z} {Surface{var2+6*22}; Layers{n_layer_z}; Recombine; }
Extrude {0,0, z} {Surface{var2+7*22}; Layers{n_layer_z}; Recombine; }

var1 = 256;
var2 = 278;
Extrude {0,0, PML_w} {Surface{var1}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{var2}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{var2+1*22}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{var2+2*22}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{var2+3*22}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{var2+4*22}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{var2+5*22}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{var2+6*22}; Layers{PML_layers}; Recombine; }
Extrude {0,0, PML_w} {Surface{var2+7*22}; Layers{PML_layers}; Recombine; }


If (0==1)
//Physical Surface ( 1 ) = { 608, 630, 476, 586, 454, 489, 520, 542, 564  };
//Physical Volume ( 1  ) = { 1,2,3,4,5,6,7,8,9  };
Physical Volume ( 1  ) = { 10 };
Physical Volume ( 2  ) = { 20 };
Physical Volume ( 3  ) = { 21 };
Physical Volume ( 4  ) = { 22 };
Physical Volume ( 5  ) = { 27 };
Physical Volume ( 6  ) = { 19 };
Physical Volume ( 7  ) = { 23 };
Physical Volume ( 8  ) = { 26 };
Physical Volume ( 9  ) = { 25 };

//Physical Volume ( 2  ) = { 10,11,12,13,14,15,16,17,18 };
Physical Volume ( 10  ) = { 24 };
Physical Volume ( 11  ) = { 11 };
Physical Volume ( 12  ) = { 12 };
Physical Volume ( 13  ) = { 13 };
Physical Volume ( 14  ) = { 18 };
Physical Volume ( 15  ) = { 14 };
Physical Volume ( 16  ) = { 17 };
Physical Volume ( 17  ) = { 16 };
Physical Volume ( 18  ) = { 15 };

Physical Volume ( 19 ) = { 2  };
Physical Volume ( 20 ) = { 3  };
Physical Volume ( 21 ) = { 4  };
Physical Volume ( 22 ) = { 9  };
Physical Volume ( 23 ) = { 1  };
Physical Volume ( 24 ) = { 5  };
Physical Volume ( 25 ) = { 8  };
Physical Volume ( 26 ) = { 7  };
Physical Volume ( 27 ) = { 6  };
EndIf


Physical Volume ( 1  ) = { 1 };
Physical Volume ( 2  ) = { 2 };
Physical Volume ( 3  ) = { 3 };
Physical Volume ( 4  ) = { 4 };
Physical Volume ( 5  ) = { 5 };
Physical Volume ( 6  ) = { 6 };
Physical Volume ( 7  ) = { 7 };
Physical Volume ( 8  ) = { 8 };
Physical Volume ( 9  ) = { 9 };

Physical Volume ( 10  ) = { 10 };
Physical Volume ( 11  ) = { 11 };
Physical Volume ( 12  ) = { 12 };
Physical Volume ( 13  ) = { 13 };
Physical Volume ( 14  ) = { 14 };
Physical Volume ( 15  ) = { 15 };
Physical Volume ( 16  ) = { 16 };
Physical Volume ( 17  ) = { 17 };
Physical Volume ( 18  ) = { 18 };

Physical Volume ( 19 ) = { 19  };
Physical Volume ( 20 ) = { 20  };
Physical Volume ( 21 ) = { 21  };
Physical Volume ( 22 ) = { 22  };
Physical Volume ( 23 ) = { 23  };
Physical Volume ( 24 ) = { 24  };
Physical Volume ( 25 ) = { 25  };
Physical Volume ( 26 ) = { 26  };
Physical Volume ( 27 ) = { 27  };


Physical Surface ( 1 ) = { 1  };
Physical Surface ( 2 ) = { 2, 67, 79  };
Physical Surface ( 3 ) = { 3,89  };
Physical Surface ( 4 ) = { 4, 111, 115  };
Physical Surface ( 5 ) = { 5, 137  };
Physical Surface ( 6 ) = { 6,159, 163  };
Physical Surface ( 7 ) = { 7, 185  };
Physical Surface ( 8 ) = { 8, 207, 211  };
Physical Surface ( 9 ) = { 9, 233  };

//Physical Surface ( 11 ) = { 277, 265 };

Physical Surface ( 11 ) = { 277, 265 };
Physical Surface ( 12 ) = { 287 };
Physical Surface ( 13 ) = { 309, 313 };
Physical Surface ( 14 ) = { 335 };
Physical Surface ( 15 ) = { 357, 361 };
Physical Surface ( 16 ) = { 383 };
Physical Surface ( 17 ) = { 405, 409 };
Physical Surface ( 18 ) = { 431 };

Physical Surface ( 19 ) = { 454 };
Physical Surface ( 20 ) = { 476, 463, 475 };
Physical Surface ( 21 ) = { 485, 498 };
Physical Surface ( 22 ) = { 507, 511, 520 };
Physical Surface ( 23 ) = { 533, 542 };
Physical Surface ( 24 ) = { 564, 559, 555 };
Physical Surface ( 25 ) = { 586, 581 };
Physical Surface ( 26 ) = { 603, 607, 608 };
Physical Surface ( 27 ) = { 629, 630 };

