// Gmsh project created on Mon Jan 13 10:20:25 2025
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
Line(203)={7,12};

//Ligne du niveu d'eau
DefineConstant [ eau = {0.25,Min 0,Max H-h,Step 0.05, Name "Parametre/niveau eau"} ] ;
Point(8)={(H-h-eau)*lb/(H-h)+l,h+eau,0};
Point(11)={(-1)*(H-h-eau)*lb/(H-h)+l+b+2*lb,h+eau,0};
Point(3)={0,h+eau,0};
Point(16)={L,h+eau,0};

Line(300)={8,11};  

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
Curve loop(1)={1,2,3,4,5,6,203,12,13,14,15,16,17,18,19};
Plane Surface(1)={1};

