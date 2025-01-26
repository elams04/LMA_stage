# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:28:17 2025

@author: Utilisateur
"""
import numpy as np
import math
import SEM_bis as sem

dx=0.05
r=2

mat=np.eye(1)
rho=np.eye(1)
Longueur=10
duree=30
L=int(Longueur/dx) #nombre d'élements mais L+1 noeuds

c=1


f=0.5
periode=1/f
a=Longueur/2    #endroit où de la source
exemplepoint=1
amplitude=10
N_point=L*r+1
instant=20

x,w= sem.lglnodes(r);
phi,dphi =sem.ShapeFunctionLagrange(x,method_type='SEM')
Mesh,X,w=sem.maillage(L,r,dx)
K_local,M_local=sem.fabricmatrice(X,w,phi,dphi,rho,mat)   
K=sem.assemble(K_local,N_point) #Matrice de rigidité
M=sem.assemble(M_local,N_point)
M_inverse=np.diag(1/np.diag(M))

def erreur(dt):
    instant_=int(instant/dt)
    ntemps=int(duree/dt) +1
    t=np.linspace(0,duree,ntemps)
    ua=sem.condition(t,amplitude,f)
    
    U=sem.FDM(duree,Longueur,r,f,a,ua,K,M_inverse,dt,dx,c)

