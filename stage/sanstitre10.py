# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 16:09:46 2024

@author: Utilisateur
"""

import numpy as np
import matplotlib.pyplot as plt

def wave_displacement(t_values, x_values, a, c, omega, T, L):
    """
    Calcule le déplacement u(x, t) pour un ensemble de points x donné en fonction du temps dans un milieu fini.

    Paramètres :
    - t_values : np.ndarray, tableau des instants de temps
    - x_values : np.ndarray, tableau des positions où on calcule le déplacement
    - a : float, position de la source
    - c : float, vitesse de propagation de l'onde
    - omega : float, fréquence angulaire de la source
    - T : float, durée pendant laquelle la source est active
    - L : float, longueur du domaine fini

    Retourne :
    - u_values : np.ndarray, matrice des valeurs de u(x, t) (dimension len(x_values) x len(t_values))
    """
    u_values = np.zeros((len(x_values), len(t_values)))
    
    for i, t in enumerate(t_values):
        for j, x in enumerate(x_values):
            # Temps retardé
            t_retard = t - abs(x - a) / c
            
            # Condition pour que le signal atteigne x et que la source soit active
            if t_retard > 0 and t_retard < T:
                u_values[j, i] = (1 / (2 * c)) * np.sin(omega * t_retard)
    
    return u_values

# Paramètres de l'équation
a = 5.0          # Position de la source
t_max = 30     # Durée totale de la simulation
c = 2.0          # Vitesse de propagation
omega = 2 * np.pi*0.5 # Fréquence angulaire de la source
T = 2          # Durée d'activation de la source
L = 10.0         # Longueur du domaine fini

# Points où on calcule le déplacement
x_values = np.linspace(0, L, 100)

# Vérification de la validité de a
if a < 0 or a > L:
    raise ValueError("La position a doit être dans le domaine [0, L].")

# Discrétisation temporelle
dt = 0.01
t_values = np.arange(0, t_max, dt)

# Calcul du déplacement
u_values = wave_displacement(t_values, x_values, a, c, omega, T, L)

# Visualisation du déplacement en un point spécifique
x_specific = 10.0
if x_specific < 0 or x_specific > L:
    raise ValueError("La position x_specific doit être dans le domaine [0, L].")

idx_specific = np.argmin(np.abs(x_values - x_specific))
u_specific = u_values[idx_specific, :]

plt.figure(figsize=(10, 6))
plt.plot(t_values, u_specific, label=f'Déplacement à x={x_specific}')
plt.axvline(x=(x_specific-a)/c, color='r', linestyle='--', label="Temps d'arrivée de l'onde")
plt.axvline(x=(x_specific-a)/c + T, color='g', linestyle='--', label="Fin de l'activation")
plt.title("Déplacement en un point en fonction du temps")
plt.xlabel("Temps t")
plt.ylabel("Déplacement u(x, t)")
plt.legend()
plt.grid()
plt.show()

# Visualisation du déplacement en fonction de l'espace et du temps (graphique type chaleur)
plt.figure(figsize=(12, 8))
plt.imshow(u_values, extent=[0, t_max, 0, L], origin='lower', aspect='auto', cmap='hot')
plt.colorbar(label="Déplacement u(x, t)")
plt.title("Carte de chaleur du déplacement en fonction de l'espace et du temps")
plt.xlabel("Temps t")
plt.ylabel("Position x")
plt.show()
