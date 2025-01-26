# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:28:17 2025

@author: Utilisateur
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons

import SEM_bis as sem
from solution_exact_bis import matrice_depl_exact
from SEM_approche_modifie_bis import approcheModifié

dx=0.1
r=2
k=6

mat=np.eye(1)
rho=np.eye(1)
Longueur=10
duree=30
L=int(Longueur/dx) #nombre d'élements mais L+1 noeuds

c=1


f=0.5
omega=np.pi*2*f
T=1/f
periode=1/f
a=Longueur/2    #endroit où de la source
exemplepoint=1
amplitude=10
N_point=L*r+1
instant=20
seuil=20 #valeur prise si valeur infini

x,w= sem.lglnodes(r);
phi,dphi =sem.ShapeFunctionLagrange(x,method_type='SEM')
Mesh,X,w=sem.maillage(L,r,dx)
K_local,M_local=sem.fabricmatrice(X,w,phi,dphi,rho,mat)   
K=sem.assemble(K_local,N_point) #Matrice de rigidité
M=sem.assemble(M_local,N_point)
M_inverse=np.diag(1/np.diag(M))


def integral(f):
    valeur=0
    c=(f[0]+f[N_point-r-1])/2
    for i in range(r,N_point-r,r):
        valeur= valeur + f[i]
    valeur=dx*(valeur + c)
    return valeur

    

    

def erreur(dt,dx):
    instant_=int(instant/dt)
    ntemps=int(duree/dt) +1
    t=np.linspace(0,duree,ntemps)
    ua=sem.condition(t,amplitude,f)
    
    U=sem.FDM(duree,Longueur,r,f,a,ua,K,M_inverse,dt,dx,c)
    U_exact=matrice_depl_exact(t, Mesh, a, c, amplitude, omega, T, Longueur)
    F=(U[:,instant_]-U_exact[:,instant_])**2
    G=(U_exact[:,instant_])**2
    erreur=integral(F)/integral(G)
    return erreur

def erreur_2(k,dt,dx):
    instant_=int(instant/dt)
    ntemps=int(duree/dt) +1
    t=np.linspace(0,duree,ntemps)
    ua=sem.condition(t,amplitude,f)
    
    U=approcheModifié(duree, Longueur, r, f, a, ua, K, M_inverse, dt, dx, c, k)
    U_exact=matrice_depl_exact(t, Mesh, a, c, amplitude, omega, T, Longueur)
    F=(U[:,instant_]-U_exact[:,instant_])**2
    G=(U_exact[:,instant_])**2
    erreur=integral(F)/integral(G)
    return erreur

# # pour avoir juste ordre 2
# temps=[0.002*i for i in range(1,22)]
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
# plt.figure(figsize=(8, 5))
# for j in range(1,k+1):
#     temps=[0.03+0.005*i for i in range(1,15)]
#     alpha=[c*dt/dx for dt in temps]
#     Erreur=[erreur_2(j,t,dx) for t in temps]
#     plt.plot(alpha,Erreur,'x',label=f'ordre {2*j}')

    
# plt.legend()
# plt.xlabel("alpha")
# plt.ylabel("Erreur")
# plt.title("Erreur en fonction du CFL FDM pour différent ordre d'approximation en temps")
# plt.grid()
# plt.show()

#Pour pouvoir enlever certaines courbes
# ordre 2 à 2k
temps=[0.03+0.001*i for i in range(0,100)]
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

