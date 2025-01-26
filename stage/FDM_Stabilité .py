# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 11:35:07 2025

@author: Utilisateur
"""
import numpy as np
from matplotlib.widgets import CheckButtons

from FDM_bis import matrice_depl_fdm
from solution_exact_bis import matrice_depl_exact
from FDM_approche_mod import matrice_depl_fdm_mod
import matplotlib.pyplot as plt
c=1
dx=0.05
k=6
Longueur=10
duree=30
exemplepoint=1
f=0.5
T=1/f
omega=2*np.pi*f
Amplitude=10 
a=Longueur/2
instant=20
seuil=20

def integral(n,f,dx):
    valeur=0
    ca=(f[0]+f[n-1])/2
    for i in range(1,n-1):
        valeur= valeur + f[i]
    valeur=dx*(valeur + ca)
    return valeur

def erreur(dt,dx):
    instant_=int(instant/dt)
    ntemps=int(duree/dt) +1
    nespace=int(Longueur/dx)+1
    t=np.linspace(0,duree,ntemps)
    X=np.linspace(0,Longueur,nespace)
    
    
    U=matrice_depl_fdm(c, duree, Longueur, a, f, Amplitude, dt, dx)
    U_exact=matrice_depl_exact(t,X,a,c,Amplitude,T,Longueur,omega)
    F=(U[:,instant_]-U_exact[:,instant_])**2
    K=(U_exact[:,instant_])**2
    erreur=integral(nespace,F,dx)/integral(nespace,K,dx)
    return erreur

def erreur_2(k,dt,dx):
    instant_=int(instant/dt)
    ntemps=int(duree/dt) +1
    nespace=int(Longueur/dx)+1
    t=np.linspace(0,duree,ntemps)
    X=np.linspace(0,Longueur,nespace)
    omega=2*np.pi*f
    
    U=matrice_depl_fdm_mod(k,c, duree, Longueur, a, f, Amplitude, dt, dx)
    U_exact=matrice_depl_exact(t,X,a,c,Amplitude,T,Longueur,omega)
    F=(U[:,instant_]-U_exact[:,instant_])**2
    K=(U_exact[:,instant_])**2
    erreur=integral(nespace,F,dx)/integral(nespace,K,dx)
    return erreur



# # pour avoir juste ordre 2
#dt varie
# temps=[0.040+0.001*i for i in range(0,30)]
# alpha=[c*dt/dx for dt in temps]
# Erreur=[]
# for i,t in enumerate (temps):
#     valeur_erreur=erreur(t,dx)
#     if valeur_erreur>=seuil:
#         abscisse=i
#         alpha=alpha[:i+1]
#         Erreur.append(seuil)
#         break
#     Erreur.append(valeur_erreur)

# dx varie
# dt=0.01
# espace=[0.005+0.001*i for i in range(0,50)]
# alpha=[c*dt/dx for dt in espace]
# Erreur=[]
# for i,x in enumerate (espace):
#     valeur_erreur=erreur(dt,x)
#     if valeur_erreur>=seuil:
#         abscisse=i
#         alpha=alpha[:i+1]
#         Erreur.append(seuil)
#     Erreur.append(valeur_erreur)
        

# # ticks=np.arange(0,alpha[abscisse],0.05)
# plt.figure(figsize=(8, 5))
# plt.plot(alpha,Erreur,'x',label='ordre 2')
# plt.xlabel("alpha")
# plt.ylabel("Erreur")
# plt.title("Erreur en fonction SEM  du CFL")
# # plt.xticks(ticks)
# plt.legend()
# plt.grid()
# plt.show()

# # ordre 2 à 2k
# temps=[0.045+0.0001*i for i in range(0,100)]
# alpha=[c*dt/dx for dt in temps]

# plt.figure(figsize=(8, 5))
# for j in range(1,k+1):
#     Erreur=[]
#     for i,t in enumerate (temps):
#         valeur_erreur=erreur_2(j,t,dx)
#         if valeur_erreur>=seuil:
#             abscisse=i
#             alpha=alpha[:i+1]
#             Erreur.append(seuil)
#             break
#         Erreur.append(valeur_erreur)
    
#     plt.plot(alpha,Erreur,'x',label=f'ordre {2*j}')

 
# plt.legend()

# plt.xlabel("alpha")
# plt.ylabel("Erreur")
# plt.title("Erreur en fonction du CFL FDM pour différent ordre d'approximation en temps")
# plt.grid()
# plt.show()

##
# ordre 2 à 2k
temps=[0.043+0.001*i for i in range(0,100)]
alpha=[c*dt/dx for dt in temps]

# Préparation des données
erreurs_par_ordre = []
labels = []
for j in range(1, k + 1):
    Erreur = []
    for t in temps:
        valeur_erreur = erreur_2(j, t,dx)
        if valeur_erreur >= seuil:
            Erreur.append(seuil)
            break
        Erreur.append(valeur_erreur)
    erreurs_par_ordre.append(Erreur)
    labels.append(f'ordre {2 * j}')

# Création de la figure
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25)  # Laisser de l'espace pour les widgets

# Tracer les courbes
courbes = []
for j, Erreur in enumerate(erreurs_par_ordre):
    curve, = ax.plot(alpha[:len(Erreur)], Erreur, 'x', label=labels[j])
    courbes.append(curve)

# Configuration des axes
ax.set_xlabel("alpha")
ax.set_ylabel("Erreur")
ax.set_title("Erreur en fonction du CFL FDM pour différents ordres")
ax.legend()
ax.grid(True)

# Création des cases à cocher
rax = plt.axes([0.01, 0.4, 0.2, 0.3])  # Position des cases à cocher
check = CheckButtons(rax, labels, [True] * len(labels))

# Fonction pour gérer la visibilité des courbes
def toggle_visibility(label):
    index = labels.index(label)
    courbes[index].set_visible(not courbes[index].get_visible())
    plt.draw()

# Connecter les cases à cocher
check.on_clicked(toggle_visibility)

# Afficher le graphique
plt.show()


    

