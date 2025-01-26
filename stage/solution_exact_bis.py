# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 09:40:58 2024

@author: Utilisateur
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation



# Fonction de déplacement u(x, t) dans un domaine fini
def displacement_finite_domain(t, x, a, c, A, omega, T, L):
    delay = np.abs(x - a) / c  # Temps de retard
    u = np.zeros_like(t)
    
    # Onde seulement si t > retard et t < retard + T, avec réflexions aux bords
    for n in range(-5, 6):  # Réflexions multiples (n = nombre de réflexions)
        x_reflected = a + (-1)**n * (x - a + n * L)  # Position réfléchie
        delay_reflected = np.abs(x_reflected - a) / c
        mask = (t >= delay_reflected) & (t <= delay_reflected + T)
        u[mask] += A * np.sin(omega * (t[mask] - delay_reflected))
    
    return u

#fonction renvoyant la matrice de deplacement (ligne=espace colone=temps)
def matrice_depl_exact(t,X,a,c,A,omega,T,L):
    U=np.zeros((len(X),len(t)))
    for j,x in enumerate (X):
        U[j,:]=displacement_finite_domain(t,x,a,c,A,omega,T,L)
    return U



if __name__ == "__main__":
    
    # Paramètres physiques
    x_0 = 1.0     # Point d'observation
    instant=20    #instant d'observation
    L = 10.0      # Longueur du domaine fini
    a = L/2       # Position de la source
    c = 1.0       # Vitesse de propagation de l'onde
    A = 10.0      # Amplitude de l'onde
    f=0.5     # frequence
    omega = 2.0 * np.pi *f  # Fréquence angulaire (ajustée)
    T = 1/f  # Durée d'émission qui corespond à une periode
    dt=0.01
    dx=0.05
    
    # Paramètres temporels et spacial
    t_max = 30.0  # Durée totale de simulation
    n_points = int(t_max/dt)+1  # Nombre de points temporels
    t = np.linspace(0, t_max, n_points) 
    X=np.linspace(0,L,int(L/dx)+1)
    
    # Calcul du déplacement au point x_0
    u = displacement_finite_domain(t, x_0, a, c, A, omega, T, L)

    # Tracer le déplacement en fonction du temps
    plt.figure(figsize=(8, 5))
    plt.plot(t, u, label=f'Déplacement en x={x_0}')
    plt.axvline(np.abs(x_0 - a) / c, color='r', linestyle='--', label='Temps de retard')
    plt.xlabel('Temps t')
    plt.ylabel('Déplacement u(x, t)')
    plt.title("Déplacement au point d'observation  dans un domaine fini")
    plt.legend()
    plt.grid()
    # # plt.savefig("slt_exact")
    plt.show()

    # # Visualisation de la solution exacte carte type chaleur
    U=matrice_depl_exact(t,X,a,c,A,omega,T,L)
    
            
            
    # # # Visualisation de la solution exacte 
    # plt.figure(figsize=(10, 6))
    # plt.title("Solution exacte de l'équation de d'Alembert")
    # plt.pcolormesh(t, X, U, shading='auto', cmap='jet')
    # plt.colorbar(label="Amplitude")
    # plt.xlabel("Temps (s)")
    # plt.ylabel("Position (x)")
    # plt.grid(True)
    # plt.show()
    
    # Tracer le déplacement en fonction de la position
    plt.figure(figsize=(8, 5))
    plt.plot(X, U[:,int(instant/dt)], label=f'Déplacement en t={x_0}')
    plt.xlabel('Temps t')
    plt.ylabel('Déplacement u(x, t)')
    plt.title("Déplacement sur tout le domaine en un temps donées")
    plt.legend()
    plt.grid()
    # # plt.savefig("slt_exact")
    plt.show()
    

    #Animation

    # fig=plt.figure()
    # plt.xlabel("longueur")
    # plt.ylabel("deplacement")
    # plt.title("probleme 1D propagation onde")
    # line,=plt.plot([],[])
    # plt.xlim(0,L)
    # plt.ylim(-A*5,A*5)


    # def animate(i):
    #      line.set_data(X,U[:,i])
    #      return line,
    # ani = animation.FuncAnimation(fig, animate, frames=range(0, n_points, 1), interval=10, blit=True, repeat=False)

    # # Enregistrement de l'animation
    # ani.save("propagation_ondes_solution_exacte_bis.mp4", writer="ffmpeg", fps=30)
    # plt.show()
