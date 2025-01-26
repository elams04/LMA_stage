# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 13:44:23 2024

@author: Utilisateur
"""

import pickle
import matplotlib.pyplot as plt

# Charger le fichier contenant la figure
with open('Stabilité_alpha.pkl', 'rb') as f:
    fig = pickle.load(f)

# Accéder aux axes du graphique
ax = fig.axes[0]  # Supposons qu'il y a un seul subplot dans la figure
ax.set_ylim(0,30)
# Lister les courbes (scatter plots stockés dans ax.collections)
collections = ax.collections

# Supprimer les 5 dernières courbes en utilisant les indices
for i in range():
    if len(collections) > 0:  # Vérifier qu'il reste des courbes
        collections[-1].remove()  # Supprimer la dernière courbe

# Afficher le graphique modifié
plt.show(fig)
