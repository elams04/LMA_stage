# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 14:37:40 2024

@author: Utilisateur
"""
import gmsh
import re

# Initialisation d'un dictionnaire pour stocker les variables
variables = {}

# Lecture du fichier .geo
with open("Barrage_V3.geo", "r") as fichier:
    for ligne in fichier:
        # Extraction des lignes de type "clé = valeur;"
        match = re.match(r"(\w+)\s*=\s*([\d.]+);", ligne)
        if match:
            clé, valeur = match.groups()
            variables[clé] = float(valeur)  # Convertit en float si c'est un nombre

# Utilisation des variables
print(variables)  # Affiche {'a': 10.0, 'b': 20.0}
print("Valeur de a :", variables["c"])

vpSol=6000
vsSol=3400
dSol=2650
vpBeton=3992
vsBeton=2444
dBeton=2370
vpEau=1400
vsEau=0
dEau=1000
point_plany=[650,282,240,514,543,264,617,217,266]

# Initialiser Gmsh
gmsh.initialize()

# Charger le fichier .geo
gmsh.open("Barrage_V3.geo")

# Synchroniser pour s'assurer que le modèle est bien chargé
gmsh.model.geo.synchronize()

# Récupérer les informations d'un point en utilisant son identifiant


coords={p: gmsh.model.getValue(0,p, []) for p in point_plany } 


# print(f"Coordonnées du point {point_id} :", coords)

# Terminer la session Gmsh
gmsh.finalize()


Materials=[["S",vpSol,vsSol,dSol],["S",vpBeton,vsBeton,dBeton],["F",vpEau,vsEau,dEau]]
for i in range(18):
    if i==12:
        Materials.append(["L",0,0,0])
    else:
        Materials.append(["P",0,0,0])
        

def material_write(nomfichier,nmat,Mater,npow,A,coords,dx,dy,dz):
    with open(nomfichier, "w+") as fichier:
        fichier.write(f"{nmat}\n") #nombre de matériaux
        #ecriture des propriétés des matériaux
        for typ,vp,vs,d in Mater:
            fichier.write(f"{typ} {vp:.8f} {vs:.8f} {d:.8f} 0.000000 0.000000\n")
        
        fichier.write("# PML properties\n")
        fichier.write("# npow,Apow,posX,widthX,posY,widthY,posZ,widthZ,mat\n")
        
        Pml=[
            [0,0,0,0,coords[650][2],dz],
            [coords[514][0],-dx,0,0,coords[514][2],dz],
            [coords[543][0],-dx,coords[543][1],-dy,coords[543][2],dz],
            [0,0,coords[617][1],-dy,coords[617][2],dz],
            [coords[617][0],dx,coords[617][1],-dy,coords[617][2],dz],
            [coords[650][0],dx,0,0,coords[650][2],dz],
            [0,0,0,0,coords[282][2],-dz],
            [coords[240][0],-dx,0,0,coords[240][2],-dz],
            [coords[217][0],dx,coords[217][1],-dy,coords[217][2],-dz],
            [0,0,coords[217][1],-dy,coords[217][2],-dz],
            [coords[264][0],dx,coords[264][1],-dy,coords[264][2],-dz],
            [coords[282][0],dx,0,0,coords[282][2],-dz],
            [0,0,0,0,coords[266][2],-dz],
            [coords[240][0],-dx,0,0,0,0],
            [coords[217][0],-dx,0,0,0,0],
            [coords[264][0],dx,0,0,0,0],
            [coords[282][0],dx,0,0,0,0],
            [0,0,coords[264][0],-dy,0,0]
            ]
        i=1
        for posX, widthX, posY, widthY, posZ, widthZ in Pml:
            if i==13:
                fichier.write(f"{npow} {A:.6f} {posX:.8f} {widthX:.8f} {posY:.8f} {widthY:.8f} {posZ:.8f} {widthZ:.8f} 3\n")
            else:
                fichier.write(f"{npow} {A:.6f} {posX:.8f} {widthX:.8f} {posY:.8f} {widthY:.8f} {posZ:.8f} {widthZ:.8f} 1\n")
            i=i+1
    return
material_write("material.input",21,Materials,2,10,coords,variables["c"],variables["cb"],variables["pml"])
        