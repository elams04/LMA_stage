# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 14:41:57 2024

@author: Utilisateur
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def solution_reflechie_temporaire(x, a, c, L, A, omega, T, t_values):
    """
    Calcule la solution de l'onde avec réflexions et source active pendant une durée limitée.
    
    Parameters:
    x (float): Position où calculer la solution.
    a (float): Position de la source Dirac.
    c (float): Vitesse de propagation de l'onde.
    L (float): Longueur du domaine.
    A (float): Amplitude du second membre.
    omega (float): Fréquence angulaire.
    T (float): Durée pendant laquelle la source est active.
    t_values (list or np.array): Liste des instants t.
    
    Returns:
    np.array: Valeurs de u(x,t) pour chaque t.
    """
    u_values = np.zeros_like(t_values)
    for n in range(-10, 10):  # Sommation sur les réflexions
        x_image = 2 * n * L + ((-1)**n) * a
        tau = t_values - np.abs(x - x_image) / c
        active_source = (tau >= 0) & (tau <= T)
        u_values += (A / (2 * c)) * np.sin(omega * tau) * active_source
    return u_values

def Matrice_deplacement(Mesh,a,c,L,A,omega,T,t_values):
    N_espace=len(Mesh)-1
    N=len(t_values)-1
    U=np.zeros((N_espace+1,N+1))
    
    for i in range(N_espace+1):
        U[i,:]=solution_reflechie_temporaire(Mesh[i],a,c,L,A,omega,T,t_values)
    return U


# Paramètres
L = 10.0        # Longueur du domaine
Duree=30     #Duree de simulation
a = L/2         # Position de la source Dirac
c = 1         # Vitesse de propagation
A = 100
f=0.5      # Amplitude
omega = 2*np.pi*f    # Fréquence angulaire
T = 1/f                 # Durée de l'activation de la source
x = 1.0       # Point d'observation
dt=0.001
dx=0.1 
N=int(Duree/dt)
N_espace=int(L/dx)     
t_values = np.linspace(0, Duree,N+1 )
Mesh=[i*dx for i in range(N_espace+1)]
# Calcul de la solution
u_values = solution_reflechie_temporaire(x, a, c, L, A, omega, T, t_values)



# U=Matrice_deplacement(Mesh,a,c,L,A,omega,T,t_values)


# Visualisation
plt.plot(t_values, u_values)
plt.title(f"Solution de l'onde en x={x} ")
plt.xlabel("Temps t")
plt.ylabel("Amplitude u(x,t)")
plt.grid(True)
plt.show()

# # Visualisation de la solution exacte
# plt.figure(figsize=(10, 6))
# plt.title("Solution différence finies de l'équation de d'Alembert")
# plt.pcolormesh(t_values, Mesh, U, shading='auto', cmap='jet')
# plt.colorbar(label="Amplitude")
# plt.xlabel("Temps (s)")
# plt.ylabel("Position (x)")
# plt.grid(True)
# plt.show()

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# x, y = np.meshgrid(t_values, Mesh)
# surf = ax.plot_surface(x, y, U, cmap='viridis', edgecolor='none')
# ax.set_title("Déplacement de l'onde dans l'espace-temps")
# ax.set_xlabel("Temps t")
# ax.set_ylabel("Position x")
# ax.set_zlabel("Amplitude u(x,t)")
# plt.show()

# fig=plt.figure()
# plt.xlabel("longueur")
# plt.ylabel("deplacement")
# plt.title("probleme 1D propagation onde")
# line,=plt.plot([],[])
# plt.xlim(0,L)
# plt.ylim(-A,A)


# def animate(i):
#      line.set_data(Mesh,U[:,i])
#      return line,
# ani = animation.FuncAnimation(fig, animate, frames=range(0, N, 1), interval=10, blit=True, repeat=False)

# # Enregistrement de l'animation
# ani.save("propagation_ondes_solution_exacte.mp4", writer="ffmpeg", fps=30)
# plt.show()