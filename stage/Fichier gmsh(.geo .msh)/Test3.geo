// Gmsh project created on Wed Oct 09 09:26:47 2024
Point(1) = {0, 0, 0};
Point(2) = {0, 1, 0};
Point(3) = {0, 1.5, 0};
Point(4) = {1, 1.5, 0};
Point(5) = {1, 1, 0};
Point(6) = {1, 0, 0};

// Création des lignes
i = 1;
For(1:5)
  Line(i) = {i, i+1};
  i++;
EndFor
Line(i) = {6, 1};
Line(7) = {2, 5};

// Création des surfaces
Curve Loop(1) = {1, 7, 5, 6};
Plane Surface(1) = {1};
Transfinite Surface {1} = {1, 2, 5, 6};
Recombine Surface {1};

Curve Loop(2) = {2, 3, 4, -7};
Plane Surface(2) = {2};
Transfinite Surface {2} = {2, 3, 4, 5};
Recombine Surface {2};

// Définition d'une boucle de surfaces
Surface Loop(1) = {1, 2};

// Extrusion combinée des surfaces via Surface Loop
Extrude {0, 0, -2} {
  Surface {1};  // Extrude les deux surfaces
  Layers{20, 1};
  Recombine;
}

Coherence;