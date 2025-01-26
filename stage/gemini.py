# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 16:53:49 2024

@author: Utilisateur
"""

import numpy as np
import matplotlib.pyplot as plt

def dAlembert_1D(L, T, dx, dt, A, omega, a):
    """
    Résout l'équation de d'Alembert en 1D avec un second membre de type Dirac.

    Args:
        L: Longueur du domaine spatial.
        T: Durée de la simulation.
        dx: Pas d'espace.
        dt: Pas de temps.
        A: Amplitude de la source.
        omega: Pulsation de la source.
        a: Position de la source.
    """

    # Vérification de la condition de stabilité de Courant-Friedrichs-Lewy (CFL)
    c = 1  # Vitesse de propagation de l'onde (à ajuster selon le problème)
    if c * dt / dx > 1:
        raise ValueError("La condition de stabilité CFL n'est pas respectée.")

    # Nombre de points d'espace et de temps
    nx = int(L/dx) + 1
    nt = int(T/dt) + 1

    # Création de la matrice de solution
    U = np.zeros((nt, nx))

    # Conditions initiales
    U[0, :] = 0  # u(x, 0) = 0
    U[1, :] = 0  # u_t(x, 0) = 0

    # Conditions aux limites
    U[:, 0] = 0  # u(0, t) = 0
    U[:, -1] = 0  # u(L, t) = 0

    # Schéma aux différences finies
    for n in range(1, nt-1):
        for i in range(1, nx-1):
            # Calcul du second membre
            f = A * np.sin(omega * n * dt) * (n * dt <= T) * (np.abs(i*dx - a) < dx/2)
            U[n+1, i] = 2 * U[n, i] - U[n-1, i] + (c * dt / dx)**2 * (U[n, i+1] - 2*U[n, i] + U[n, i-1]) + f * dt**2

    return U

# # Paramètres
# L = 10
# T = 30
# dx = 0.1
# dt = 0.1
# A = 100
# omega = 2*np.pi*0.5
# a = 5

# # Résolution
# U = dAlembert_1D(L, T, dx, dt, A, omega, a)

# # Visualisation
# x = np.linspace(0, L, U.shape[1])
# t = np.linspace(0, T, U.shape[0])
# X, T = np.meshgrid(x, t)

# # plt.figure(figsize=(10, 6))
# # plt.pcolormesh(X, T, U)
# # plt.xlabel('x')
# # plt.ylabel('t')
# # plt.title("Solution de l'équation de d'Alembert")
# # plt.colorbar()
# # plt.show()

# plt.figure()

# plt.plot(t,U[:,int(1//dx)],"b",label="amplitude en x=")
# plt.xlabel("temps")
# plt.ylabel("déplacement")
# plt.legend()
# plt.grid(True)
# plt.show()