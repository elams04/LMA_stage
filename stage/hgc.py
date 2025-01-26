# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 09:03:11 2024

@author: Utilisateur
"""

import pickle
import matplotlib.pyplot as plt

# Charger la figure depuis le fichier pickle
with open("Stabilité_alpha.pkl", "rb") as f:
    fig = pickle.load(f)

# Récupérer l'axe principal
ax = fig.axes[0]

# Liste des ordres à afficher
ordre_a_afficher = [2, 4, 6]

# Filtrer et masquer les courbes
for line in ax.get_lines():
    label = line.get_label()
    if "ordre=" in label:
        ordre = int(label.split("=")[1])
        if ordre not in ordre_a_afficher:
            line.remove()  # Supprimer la courbe si elle ne correspond pas

# Mettre à jour la légende
ax.legend(loc='best')  # Générer une nouvelle légende avec les courbes restantes

# Réafficher la figure modifiée
plt.show()

# Optionnel : sauvegarder la figure modifiée
fig.savefig("Stabilité_alpha_modifiée.png", dpi=300)