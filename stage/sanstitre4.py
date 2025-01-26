# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 10:04:58 2024

@author: Utilisateur
"""

import matplotlib.pyplot as plt

import pickle

# Charger la figure
with open('Stabilité_alpha.pkl', 'rb') as f:
    fig = pickle.load(f)

# Afficher ou modifier la figure
ax = fig.gca()

lines = ax.get_lines()  # Liste des courbes (objets Line2D)
print(lines)  # Affiche les courbes présentes  # Modifier les limites
for line in lines:
    print(line.get_label())  # Affiche le label de la courbe
    print(line.get_color()) 
