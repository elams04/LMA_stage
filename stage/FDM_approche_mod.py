# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 10:16:02 2025

@author: Utilisateur
"""
import numpy as np
import math
import matplotlib.pyplot as plt


def condition(t,A,f):
    omega=f*2*np.pi
    T=1/f
    valeur=A*np.sin(omega*t) *(t<=T)
    return valeur



def matrice_depl_fdm_mod(k,c,duree,L,a,f,Amplitude,dt,dx):
    alpha=c*dt/dx
    lvecteur=int(L/dx)
    cvecteur=int(duree/dt)
    a_=int(a/dx)+1
    gamma=alpha**2
    T=1/f
    ua=[condition(i*dt,Amplitude,f) for i in range(int(duree/dt)+1)]
    U=np.zeros((lvecteur+1,cvecteur+1))
    #Pour k=1 fdm classique 
    if k>1:
        coef_correction_temporel=[] #Liste contenant les coef correctif de la dérvée temporelle
        coef_derive_spatial=[] 
        for j in range (2,k+1):
            coef_correction_temporel.append((dt*c)**(2*j)/math.factorial(2*j)) 
            a=np.zeros((2*j+1,2*j+1))
            b=np.zeros((2*j+1,1))
            b[2*j,0]=1
            num_colone=0
            for m in range (-j,j+1):        
                colone=[m**e for e in range(0,2*j+1)]
                a[:,num_colone]=colone
                num_colone=num_colone+1
            coef_derive_spatial.append((math.factorial(2*j))*np.linalg.solve(a,b))
           #coefs de dérivé spatial d'ordre 2j obtenu en resolvant ax=b stocker dans liste
           
       #Resolution du prblm
        for n in range(2,cvecteur+1):
            for j in range(1,lvecteur):
                U[j,n]=2*U[j,n-1]-U[j,n-2] + gamma*(U[j+1,n-1]-2*U[j,n-1]+U[j-1,n-1])
                for l in range(k-1):
                    if j<lvecteur-l-1 and j>=l+2:
                        U[j,n] = U[j,n]+(coef_correction_temporel[l]*np.array(coef_derive_spatial[l]).reshape(1,2*(l+2)+1) @ U[j-l-2:j+l+3,n]).item()
                
                if n*dt<=T and j==a_:
                    U[a_,n]=ua[n]
                    
            U[0, n] = U[1, n]
            U[lvecteur, n] = U[lvecteur-1, n]
    else:
        for n in range(2,cvecteur+1):
            for j in range(1,lvecteur):
                U[j,n]=2*U[j,n-1]-U[j,n-2] + gamma*(U[j+1,n-1]-2*U[j,n-1]+U[j-1,n-1])
                if n*dt<=T and j==a_:
                    U[a_,n]=ua[n]
                    
            U[0, n] = U[1, n]
            U[lvecteur, n] = U[lvecteur-1, n]
        
    return U


if __name__ == "__main__":
    c=1
    dt=0.01
    dx=0.05
    k=1
    alpha=c*dt/dx
    Longueur=10
    duree=30
    exemplepoint=1
    f=0.5
    Amplitude=20
    a=Longueur/2 #emplacement de la source
    U=matrice_depl_fdm_mod(k,c,duree,Longueur,a,f,Amplitude,dt,dx)
    t=np.linspace(0,duree,int(duree/dt)+1)
    X=np.linspace(0,Longueur,int(Longueur/dx)+1)
    

    # plt.plot(t,U[int(exemplepoint/dx),:],"b",label=f'amplitude en x={exemplepoint}')
    # plt.xlabel("temps")
    # plt.ylabel("déplacement")
    # plt.title(f'Solution FDM ordre={2*k} en x={exemplepoint} ')
    # plt.legend()
    # plt.grid(True)
    # plt.show()

    # # Visualisation de la solution exacte carte type chaleur
    # plt.figure(figsize=(10, 6))
    # plt.title("Solution exacte de l'équation de d'Alembert")
    # plt.pcolormesh(t, X, U, shading='auto', cmap='jet')
    # plt.colorbar(label="Amplitude")
    # plt.xlabel("Temps (s)")
    # plt.ylabel("Position (x)")
    # plt.grid(True)
    # plt.show()

    #Visualisation vidéo

    # fig=plt.figure()
    # plt.xlabel("longueur")
    # plt.ylabel("deplacement")
    # plt.title("probleme 1D propagation onde")
    # line,=plt.plot([],[])
    # plt.xlim(0,Longueur)
    # plt.ylim(-Amplitude*5,Amplitude*5)


    # def animate(i):
    #      line.set_data(X,U[:,i])
    #      return line,
    # ani = animation.FuncAnimation(fig, animate, frames=range(0, int(duree/dt)+1, 1), interval=10, blit=True, repeat=False)

    # # Enregistrement de l'animation
    # ani.save("propagation_ondes_fdm_bis.mp4", writer="ffmpeg", fps=30)
    # plt.show()

