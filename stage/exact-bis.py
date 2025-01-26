# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 11:33:26 2024

@author: Utilisateur
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# Paramètres
c = 1
dt = 0.1
dx = 0.1
Longueur = 10
duree = 30
f = 0.5
Amplitude = 10
a = Longueur / 2
x = 1

# Fonction exacte
def u_exact(x, t, a, c, A, f):
    omega = f * np.pi * 2
    u = A / (2 * c) * np.sin(omega * t) * np.cos(omega * abs(x - a) / c)
    return u

def U_exact(T, X, a, c, A, f):
    U = np.zeros((len(X), len(T)))
    for i, t in enumerate(T):
        for j, x in enumerate(X):
            U[j, i] = u_exact(x, t, a, c, A, f)
    return U

# Grilles
temps = np.linspace(0, duree, int(duree / dt) + 1)
X = np.linspace(0, Longueur, int(Longueur / dx) + 1)

# Calcul des solutions
u = [u_exact(x, t, a, c, Amplitude, f) for t in temps]
U = U_exact(temps, X, a, c, Amplitude, f)

# Visualisation de la solution pour une position donnée
plt.figure()
plt.plot(temps, u)
plt.title("Solution exacte en x = 1")
plt.xlabel("Temps (s)")
plt.ylabel("Amplitude")
plt.grid()
plt.show()

# # Visualisation de la solution exacte en 2D
# plt.figure(figsize=(10, 6))
# plt.title("Solution exacte de l'équation de d'Alembert")
# plt.pcolormesh(temps, X, U, shading='auto', cmap='jet')  # Utilisation de shading='auto'
# plt.colorbar(label="Amplitude")
# plt.xlabel("Temps (s)")
# plt.ylabel("Position (x)")
# plt.grid(True)
# plt.show()

# fig=plt.figure()
# plt.xlabel("longueur")
# plt.ylabel("deplacement")
# plt.title("probleme 1D propagation onde")
# line,=plt.plot([],[])
# plt.xlim(0,Longueur)
# plt.ylim(-Amplitude,Amplitude)


# def animate(i):
#      line.set_data(X,U[:,i])
#      return line,
# ani = animation.FuncAnimation(fig, animate, frames=range(0, int(duree/dt)+1, 1), interval=10, blit=True, repeat=False)

# # Enregistrement de l'animation
# ani.save("propagation_ondes_exact_bis.mp4", writer="ffmpeg", fps=30)
# plt.show()