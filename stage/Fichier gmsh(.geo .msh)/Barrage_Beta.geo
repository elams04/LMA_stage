// Gmsh project created on Wed Oct 09 13:45:22 2024
// Definition des variables et du parametre du niveu de l eau
H=250; h=50; b=400; lb=25; l=70; L=2*l+2*lb+b;niveaubarrage=130;
c=12000;cb=24750;
p=10000; couchep=21;
pbarrage=15;couchepbarrage=2;
avant=15000;coucheavant=9;

a1=2; a2=10; a3=2; a4=25;
b1=5; b2=2; b3=31; b4=2; b5=5;

DefineConstant[eau={75,Min 0, Max H-h,Step 10, Name "PARAMETRE/niveau eau"} ];

//Lignes du niveau barrage
Point(201)={-c,h+niveaubarrage,0};
Point(4)={0,h+niveaubarrage,0};
Point(7)={(H-h-niveaubarrage)*lb/(H-h)+l,h+niveaubarrage,0};
Point(12)={(-1)*(H-h-niveaubarrage)*lb/(H-h)+l+b+2*lb,h+niveaubarrage,0};
Point(15)={L,h+niveaubarrage,0};
Point(202)={L+c,h+niveaubarrage,0};
i=1;
Line(201)={201,4}; Line(202)={4,7}; Line(203)={7,12}; Line(204)={12,15}; Line(205)={15,202};

//Ligne du niveu d'eau
DefineConstant [ eau = {0.25,Min 0,Max H-h,Step 0.05, Name "Parametre/niveau eau"} ] ;
Point(8)={(H-h-eau)*lb/(H-h)+l,h+eau,0};
Point(11)={(-1)*(H-h-eau)*lb/(H-h)+l+b+2*lb,h+eau,0};
Point(3)={0,h+eau,0};
Point(16)={L,h+eau,0};

Line(300)={8,11}; Line(301)={3,8}; Line(302)={11,16};

//Forme barrage + sol

Point(1)={0,0,0};
table[]=Point{1};
Printf("%lf%lf%lf",table[0],table[1],table[2]);
Point(2)={0,h,0};
Point(5)={0,H,0};
Point(6)={l,H,0};
Point(9)={l+lb,h,0};
Point(10)={l+lb+b,h,0};
Point(13)={L-l,H,0};
Point(14)={L,H,0};
Point(17)={L,h,0};
Point(18)={L,0,0};
Point(19)={l+lb+b,0,0};
Point(20)={l+lb,0,0};
i=1;
For(1:19)
   Line(i)={i,i+1};
   i++;
EndFor
Line(20)={20,1};
Line(21)={9,2}; Line(22)={9,20};
Line(23)={19,10}; Line(24)={10,17};

//Zone autour du barrage cotes
g[]=Translate {-c, 0, 0} { Duplicata{ Point{1,2,3,5}; } };
d[]=Translate {c, 0, 0} { Duplicata{ Point{18,17,16,14}; } };
Line(25)={1,g[0]}; Line(26)={2,g[1]}; Line(27)={3,g[2]}; Line(28)={5,g[3]};
Line(29)={18,d[0]}; Line(30)={17,d[1]}; Line(31)={16,d[2]}; Line(32)={14,d[3]};
Line(33)={g[3],201}; Line(34)={201,g[2]}; Line(35)={g[2],g[1]}; Line(36)={g[1],g[0]};
Line(37)={d[3],202}; Line(38)={202,d[2]}; Line(39)={d[2],d[1]}; Line(40)={d[1],d[0]};

//Bas de la zonne
b[]=Translate{0,-cb,0} {Duplicata{ Point{g[0],1,20,19,18,d[0]}; }};
Line(41)={g[0],b[0]};
Line(42)={1,b[1]};
Line(43)={20,b[2]};
Line(44)={19,b[3]};
Line(45)={18,b[4]};
Line(46)={d[0],b[5]};
i=0;
For(0:4)
  Line(i+47)={b[i],b[i+1]};
  i++;
EndFor
// definition des surfaces

Macro Distance
  p1[]=Point{t1};
  p2[]=Point{t2};
  p3[]=Point{t3};
  distance12=Sqrt((p2[0] - p1[0])^2 + (p2[1] - p1[1])^2 + (p2[2] - p1[2])^2);
  distance23=Sqrt((p2[0] - p3[0])^2 + (p2[1] - p3[1])^2 + (p2[2] - p3[2])^2);
  distanceT=Sqrt((p1[0] - p3[0])^2 + (p1[1] - p3[1])^2 + (p1[2] - p1[2])^2);
  coef[]={distance12/distanceT,distance23/distanceT};
 Return
// Surface 1-9 crevasse
mur[]={1,2,9,20 };
Curve Loop(1)={1,-21,22,20};
Transfinite Line{1,22}=a3;
Transfinite Line{21,20}=b2;
Plane Surface(1)={1};
Transfinite Surface{1}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{1};

mur[]={2,3,8,9};
Curve Loop(2)={2,301,8,21};
t1=7; t2=8; t3=9;
Call Distance;

Transfinite Line{2,8}=Floor(a2*coef[1])+1;
Transfinite Line{301}=b2;
Plane Surface(2)={2};
Transfinite Surface{2}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{2};

mur[]={3,8,7,4};
Curve Loop(3)={3,202,7,-301};
Transfinite Line{202}=b2;
Transfinite Line{3,7}=Floor(a2*coef[0]);
Plane Surface(3)={3};
Transfinite Surface{3}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{3};

mur[]={4,5,6,7};
Curve Loop(4)={4,5,6,-202};
Transfinite Line{6,4}=a1;
Transfinite Line{5}=b2;
Plane Surface(4)={4};
Transfinite Surface{4}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{4};

mur[]={19,20,9,10};
Curve Loop(5)={19,-22,9,-23};
Transfinite Line{23}=a3;
Transfinite Line{19,9}=b3;
Plane Surface(5)={5};
Transfinite Surface{5}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{5};

mur[]={10,17,18,19};
Curve Loop(6)={23,24,17,18};
Transfinite Line{17}=a3;
Transfinite Line{24,18}=b4;
Plane Surface(6)={6};
Transfinite Surface{6}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{6};

t1=10;t2=11;t3=12;
Call Distance;
mur[]={10,11,16,17};
Curve Loop(7)={10,302,16,-24};
Transfinite Line{302}=b4;
Transfinite Line{16,10}=Floor(coef[0]*a2)+1;
Plane Surface(7)={7};
Transfinite Surface{7}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{7};

mur[]={11,12,15,16};
Curve Loop(8)={11,204,15,-302};
Transfinite Line{204}=b4;
Transfinite Line{11,15}=Floor(coef[1]*a2);
Plane Surface(8)={8};
Transfinite Surface{8}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{8};

mur[]={12,13,14,15};
Curve Loop(9)={12,13,14,-204};
Transfinite Line{12,14}=a1;
Transfinite Line{13}=b4;
Plane Surface(9)={9};
Transfinite Surface{9}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{9};

// Surface 10-22 contour
mur[]={15,14,d[3],202};
Curve Loop(10)={32,37,-205,-14};
Transfinite Line{37}=a1;
Transfinite Line{32,205}=b5;
Plane Surface(10)={10};
Transfinite Surface{10}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{10};

mur[]={15,202,d[2],16};
Curve Loop(11)={205,38,-31,-15};
Transfinite Line{38}=Floor(coef[1]*a2);
Transfinite Line{31}=b5;
Plane Surface(11)={11};
Transfinite Surface{11}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{11};

mur[]={17,16,d[2],d[1]};
Curve Loop(12)={16,30,-39,-31};
Transfinite Line{39}=Floor(coef[0]*a2)+1;
Transfinite Line{30}=b5;
Plane Surface(12)={12};
Transfinite Surface{12}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{12};

mur[]={d[1],d[0],18,17};
Curve Loop(13)={30,40,-29,-17};
Transfinite Line{40}=a3;
Transfinite Line{29}=b5;
Plane Surface(13)={13};
Transfinite Surface{13}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{13};

mur[]={18,d[0],b[5],b[4]};
Curve Loop(14)={29,46,-51,-45};
Transfinite Line{45,46}=a4;
Transfinite Line{51}=b5;
Plane Surface(14)={14};
Transfinite Surface{14}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{14};

mur[]={18,19,b[3],b[4]};
Curve Loop(15)={18,44,50,-45};
Transfinite Line{44}=a4;
Transfinite Line{50}=b4;
Plane Surface(15)={15};
Transfinite Surface{15}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{15};

mur[]={19,20,b[2],b[3]};
Curve Loop(16)={19,43,49,-44};
Transfinite Line{43}=a4;
Transfinite Line{49}=b3;
Plane Surface(16)={16};
Transfinite Surface{16}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{16};


mur[]={20,1,b[1],b[2]};
Curve Loop(17)={20,42,48,-43};
Transfinite Line{42}=a4;
Transfinite Line{48}=b2;
Plane Surface(17)={17};
Transfinite Surface{17}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{17};

mur[]={g[0],b[0],b[1],1};
Curve Loop(18)={25,41,47,-42};
Transfinite Line{41}=a4;
Transfinite Line{47,25}=b1;
Plane Surface(18)={18};
Transfinite Surface{18}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{18};

mur[]={1,2,g[1],g[0]};
Curve Loop(19)={1,26,36,-25};
Transfinite Line{36}=a3;
Transfinite Line{26}=b1;
Plane Surface(19)={19};
Transfinite Surface{19}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{19};

mur[]={2,3,g[2],g[1]};
Curve Loop(20)={2,27,35,-26};
Transfinite Line{35}=Floor(a2*coef[0])+1;
Transfinite Line{27}=b1;
Plane Surface(20)={20};
Transfinite Surface{20}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{20};


mur[]={3,4,201,g[2]};
Curve Loop(21)={3,-201,34,-27};
Transfinite Line{34}=Floor(a2*coef[1]);
Transfinite Line{201}=b1;
Plane Surface(21)={21};
Transfinite Surface{21}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{21};

mur[]={4,5,g[3],201};
Curve Loop(22)={4,28,33,201};
Transfinite Line{33}=a1;
Transfinite Line{28}=b1;
Plane Surface(22)={22};
Transfinite Surface{22}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{22};

//Surface 23-24 barrage 24eau
mur[]={7,8,11,12};
Curve Loop(23)={7,300,11,-203};
Transfinite Line{300,203}=b3;
Plane Surface(23)={23};
Transfinite Surface{23}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{23};

mur[]={8,9,10,11};
Curve Loop(24)={8,9,10,-300};
Plane Surface(24)={24};
Transfinite Surface{24}={mur[0],mur[1],mur[2],mur[3]};
Recombine Surface{24};

// Extrusion

//Extrusion arrière sol
topo_arriere[]=Extrude{0,0,-p}{
 Surface{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24}; Layers{couchep,1}; Recombine;
};

//Extrusion barrage

barrage[]=Extrude{0,0,pbarrage}{
 Surface{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}; Layers{couchepbarrage,1}; Recombine;
};

//Extrusion sans barage


topo_avant[]=Extrude {0, 0, avant} {
  Surface{1292,1270, 1248, 1226, 1204, 1182 ,830, 852,874, 896, 918, 1160,1138, 940, 962, 984, 1006, 1028, 1050,1072, 1094,1116}; Layers {coucheavant,1}; Recombine;
};
Physical Volume("Sol",1)={1};
i=2;
For(2:22)
  Physical Volume("Sol")+={i};
  i++;
EndFor
Physical Volume("Eau",3)={23};
i=24;
For(24:45)
  Physical Volume("Sol")+={i};
  i++;
EndFor
Physical Volume("Barrage",2)={46,47};
i=48;
For(48:69)
  Physical Volume("Sol")+={i};
  i++;
EndFor


Color Yellow {Physical Volume{1}; }
Color Grey{Physical Volume{2}; }
Color Blue {Physical Volume{3}; }

Mesh 3;
Save "Barrage_Beta.msh";
