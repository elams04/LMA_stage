# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 14:10:05 2024

@author: Utilisateur
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D


def condition(t,A,f):
    omega=f*2*np.pi
    T=1/f
    valeur=A*np.sin(omega*t) *(t<=T)
    return valeur



def matrice_depl_fdm(c,duree,L,a,f,Amplitude,dt,dx):
    alpha=c*dt/dx
    lvecteur=int(L/dx)
    cvecteur=int(duree/dt)
    a_=int(a/dx)+1
    gamma=alpha**2
    T=1/f
    ua=[condition(i*dt,Amplitude,f) for i in range(int(duree/dt)+1)]
    U=np.zeros((lvecteur+1,cvecteur+1))
    
   
    
    #Resolution du prblm
    for n in range(2,cvecteur+1):
        if n*dt<=T:
            for j in range(1,lvecteur):
                U[j,n]=2*U[j,n-1]-U[j,n-2] + gamma*(U[j+1,n-1]-2*U[j,n-1]+U[j-1,n-1])
                if j==a_:
                    U[a_,n]=ua[n]
                
        else:
            for j in range (1,lvecteur):
                U[j,n]=2*U[j,n-1]-U[j,n-2] + gamma*(U[j+1,n-1]-2*U[j,n-1]+U[j-1,n-1])
                
        U[0, n] = U[1, n]
        U[lvecteur, n] = U[lvecteur-1, n]
        
    return U


if __name__ == "__main__":
    nombre_point_mini=200
    c=1
    dt=0.001
    dx=0.01
    alpha=c*dt/dx
    Longueur=10
    duree=30
    if int(Longueur/dx)+1<nombre_point_mini:
        print("pas assez de point en espace")
    exemplepoint=1
    instant=20
    f=0.5
    Amplitude=10 
    a=Longueur/2 #emplacement de la source
    U=matrice_depl_fdm(c,duree,Longueur,a,f,Amplitude,dt,dx)
    t=np.linspace(0,duree,int(duree/dt)+1)
    x=np.linspace(0,Longueur,int(Longueur/dx)+1)
    PA=U[int(exemplepoint/dx),:]
    
    plt.figure(figsize=(8, 5))
    plt.plot(t,U[int(exemplepoint/dx),:],"b",label=f'amplitude en x={exemplepoint}')
    plt.xlabel("temps en s")
    plt.ylabel("déplacement en m")
    plt.title(f'Solution FDM en x={exemplepoint} ')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    # # X,T=np.meshgrid(t,x)
    # # fig=plt.figure()
    # # ax=fig.add_subplot(111,projection="3d")
    # # ax.plot_surface(X, T, U, cmap='viridis', edgecolor='none')
    # # plt.show()
    
    # # Tracer le déplacement en fonction de la position
    # plt.figure(figsize=(8, 5))
    # plt.plot(x, U[:,int(instant/dt)], label=f'Déplacement en t={instant}')
    # plt.xlabel('position x en m')
    # plt.ylabel('Déplacement u(x, t)')
    # plt.title("Déplacement sur tout le domaine en un temps donées")
    # plt.legend()
    # plt.grid()
    # # # plt.savefig("slt_exact")
    # plt.show()

    # # Visualisation de la solution exacte carte type chaleur
    # plt.figure(figsize=(10, 6))
    # plt.title("Solution FDM de l'équation de d'Alembert")
    # plt.pcolormesh(t, x, U, shading='auto', cmap='jet')
    # plt.colorbar(label="Amplitude")
    # plt.xlabel("Temps (s)")
    # plt.ylabel("Position (x)")
    # plt.grid(True)
    # plt.show()

    # Visualisation vidéo

    fig=plt.figure()
    plt.xlabel("longueur")
    plt.ylabel("deplacement")
    plt.title("probleme 1D propagation onde")
    line,=plt.plot([],[])
    plt.xlim(0,Longueur)
    plt.ylim(-Amplitude*5,Amplitude*5)


    def animate(i):
         line.set_data(x,U[:,i])
         return line,
    ani = animation.FuncAnimation(fig, animate, frames=range(0, int(duree/dt)+1, 50), interval=0.01, blit=True, repeat=False)

    # Enregistrement de l'animation
    ani.save("propagation_ondes_fdm_bis.mp4", writer="ffmpeg", fps=30)
    plt.show()