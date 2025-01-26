# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:34:04 2024

@author: Utilisateur
"""
import numpy as np
import matplotlib.pyplot as plt
import math

import solution_exacte as exact
import FDM
# import SEM
# Paramètres
L = 30.0        # Longueur du domaine
Duree=60     #Duree de simulation
a = L/2         # Position de la source Dirac
c = 1         # Vitesse de propagation
A = 100
f=0.5                 # Amplitude
omega= 3.14 # Fréquence angulaire
T = 1/f                 # Durée de l'activation de la source
x = 1.0       # Point d'observation
dt=0.1
dx=0.1 
N=int(Duree/dt)
N_espace=int(L/dx)     
t_values = np.linspace(0, Duree,N+1 )
Mesh=[i*dx for i in range(N_espace+1)]
# Calcul de

def Force(a,L,duree,f,amplitude,dt,dx,t):
    lvecteur=int(L/dx)
    cvecteur=int(duree/dt)
    periode=1/f
    diracIndice=int(a/dx)
    # etablissement de la matrice contenant le sinus porte
    dirac=np.array([0 for i in range(lvecteur+1)])
    dirac[diracIndice]=1
    dirac=dirac.reshape(lvecteur+1,1)
    
    F=amplitude*np.sin(2* np.pi * t / periode) * (t <= periode)
    F=F.reshape(1,cvecteur+1)
    DiracF=np.dot(dirac,F)
    
    return DiracF

def erreurFDM(c,duree,L,f,A,dt,dx,F,t):
    omega=2*np.pi*f
    T=1/f
    U_exact=exact.Matrice_deplacement(Mesh,a,c,L,A,omega,T,t)
    U=FDM.Matrice_deplacement(c,duree,L,F,dt,dx)
    return math.log(np.linalg.norm(U_exact-U,ord=2),10)


    
    
    
Temps=[]
Erreur=[]    
for i in range (0,60):
    dt=0.05 +i*0.001
    Temps.append(dt*c/dx)
    
    N=int(Duree/dt)
    t_values = np.linspace(0, Duree,N+1 )
    F=Force(a,L,Duree,f,A*15.5,dt,dx,t_values)
    Erreur.append(erreurFDM(c,Duree,L,f,A,dt,dx,F,t_values))

plt.figure()
plt.plot(Temps,Erreur,'x',label='erreur FDM')
plt.title("Stabilité pour FDM")
plt.xlabel("alpha (CFL)")
plt.ylabel("log(norme(Erreur))")  
plt.title  

