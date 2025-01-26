# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 17:01:22 2024

@author: Utilisateur
"""

import numpy as np
import matplotlib.pyplot as plt

def dirac_approx(x, a, dx):
    """Approximation discrète de la fonction delta de Dirac centrée en x = a."""
    return np.where(np.abs(x - a) < dx / 2, 1 / dx, 0)

def solve_dalembert_1d(L, T, dx, dt, a, A, omega):
    """
    Résout l'équation de d'Alembert en 1D avec un terme source ponctuel de type Dirac.
    
    Paramètres :
        L (float) : Longueur de la corde (domaine spatial).
        T (float) : Durée de simulation (domaine temporel).
        dx (float) : Pas d'espace.
        dt (float) : Pas de temps.
        a (float) : Position du terme source.
        A (float) : Amplitude du terme source.
        omega (float) : Pulsation du terme source.
        
    Retourne :
        U (ndarray) : Matrice des déplacements U(x,t).
        x (ndarray) : Discrétisation spatiale.
        t (ndarray) : Discrétisation temporelle.
    """
    c = 1  # Vitesse de propagation (hypothèse : c = 1 pour simplifier)
    Nx = int(L / dx) + 1
    Nt = int(T / dt) + 1
    x = np.linspace(0, L, Nx)
    t = np.linspace(0, T, Nt)
    
    # Matrice des déplacements
    U = np.zeros((Nx, Nt))
    
    # Calcul du terme source (appliqué uniquement pendant une période)
    source = dirac_approx(x, a, dx) * A
    
    for n in range(1, Nt-1):
        for i in range(1, Nx-1):
            source_term = source[i] * np.sin(omega * t[n]) * (t[n] <= 2 * np.pi / omega)  # H(T-t)
            U[i, n+1] = (2 * (1 - (c * dt / dx) ** 2) * U[i, n]
                         - U[i, n-1]
                         + (c * dt / dx) ** 2 * (U[i+1, n] + U[i-1, n])
                         + dt**2 * source_term)
    
    return U, x, t

# Paramètres
L = 10        # Longueur de la corde
T = 30        # Durée de simulation
dx = 0.1      # Pas d'espace
dt = 0.1     # Pas de temps
a = L/2
f=1         # Position du Dirac
A = 100      # Amplitude du signal
omega = 2*np.pi*f # Pulsation (fréquence)

# Résolution
U, x, t = solve_dalembert_1d(L, T, dx, dt, a, A, omega)

# Visualisation
plt.imshow(U, extent=[0, T, 0, L], aspect='auto', origin='lower', cmap='seismic')
plt.colorbar(label='Amplitude')
plt.xlabel('Temps (s)')
plt.ylabel('Position (m)')
plt.title("Propagation de l'onde avec un terme source Dirac")
plt.show()

plt.figure()

plt.plot(t,U[int(1//dx),:],"b",label="amplitude en x=")
plt.xlabel("temps")
plt.ylabel("déplacement")
plt.legend()
plt.grid(True)
plt.show()