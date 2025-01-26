// Gmsh project created on Tue Oct 08 10:27:35 2024
lc=0.01;
x=0; y=0; z=0;
L=0.4;
l=0.2;
h=0.5;
Point(1)={x,y,z,lc};
Point(2)={x+l,y,z,lc};
Point(3)={x+l,y+L,z,lc};
Point(4)={x,y+L,z,lc};
Point(5)={x,y,z+h,lc};
Point(6)={x+l,y,z+h,lc};
Point(7)={x+l,y+L,z+h,lc};
Point(8)={x,y+L,z+h,lc};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(9)={8,5};
Curve Loop(2)={5,6,7,9};
Plane Surface(2)={2};
Curve Loop(1)={1,2,3,4};
Plane Surface(1)={1};
T=1;
For(1:4)
    Line(T+9)={T,T+4};
    T++;
EndFor
Curve Loop(3)={4,10,-9,-13};
Plane Surface(3)={3};
Curve Loop(4)={11,6,-12,-2};
Plane Surface(4)={4};
Curve Loop(5)={10,5,-11,-1};
Plane Surface(5)={5};
Curve Loop(6)={13,-7,-12,3};
Plane Surface(6)={6};
Surface Loop(1)={1,4,3,5,2,6};
Volume(1)={1};
Physical Surface("face pave")={1,2,3,4,5,6};
T=1;
For(1:6)
    Transfinite Curve{T}=10;
    Recombine Surface {T};
    T++;
EndFor
Mesh 3;
Mesh.Smoothing=20;