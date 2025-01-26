# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 13:34:21 2024

@author: Utilisateur
"""

import numpy as np


import matplotlib.pyplot as plt


def lglnodes(N):
    N1 = N + 1
    x = np.cos(np.pi * np.arange(N+1) / N)  # Initial guess: Chebyshev-Gauss-Lobatto nodes
    P = np.zeros((N1, N1))
    xold = 2 * np.ones_like(x)
    
    while np.max(np.abs(x - xold)) > np.finfo(float).eps:
        xold = x.copy()
        P[:, 0] = 1
        P[:, 1] = x
        for k in range(2, N1):
            P[:, k] = ((2 * k - 1) * x * P[:, k-1] - (k - 1) * P[:, k-2]) / k
        x = xold - (x * P[:, N] - P[:, N-1]) / (N1 * P[:, N])
    
    w = np.flip(2 / (N * N1 * P[:, N]**2))
    return np.flip(x), w


def ShapeFunctionLagrange(x, method_type):
    if method_type == 'FEM':
        support = np.linspace(-1, 1, len(x))
    elif method_type == 'SEM':
        support = x
    phi = LagrangePolynomial(x, support)
    dphi = DifferentialLagrangePolynomial(x, support)
    return phi, dphi

def LagrangePolynomial(evaluation_points, interpolation_points):
    evaluation_points = np.asarray(evaluation_points).reshape(-1)
    interpolation_points = np.asarray(interpolation_points).reshape(-1)
    n, m = len(evaluation_points), len(interpolation_points)
    LagrangePolyEval = np.ones((n, m))
    PolynomialConstants = np.ones(m)
    
    for k in range(m):
        for i in range(m):
            if i != k:
                PolynomialConstants[k] /= (interpolation_points[k] - interpolation_points[i])
    
    for k in range(m):
        for i in range(m):
            if i != k:
                LagrangePolyEval[:, k] *= (evaluation_points - interpolation_points[i])
        LagrangePolyEval[:, k] *= PolynomialConstants[k]
    
    return LagrangePolyEval

def DifferentialLagrangePolynomial(evaluation_points, interpolation_points):
    evaluation_points = np.asarray(evaluation_points).reshape(-1)
    interpolation_points = np.asarray(interpolation_points).reshape(-1)
    n, m = len(evaluation_points), len(interpolation_points)
    LagrangePolyEval = np.zeros((n, m))
    PolynomialConstants = np.ones(m)
    
    for k in range(m):
        for i in range(m):
            if i != k:
                PolynomialConstants[k] /= (interpolation_points[k] - interpolation_points[i])
    
    for k in range(m):
        for j in range(m):
            temp = np.ones(n)
            for i in range(m):
                if i != j and i != k:
                    temp *= (evaluation_points - interpolation_points[i])
            if j != k:
                LagrangePolyEval[:, k] += temp * PolynomialConstants[k]
    
    return LagrangePolyEval



Order = 3

x,w= lglnodes(Order);
phi,dphi =ShapeFunctionLagrange(x,method_type='SEM')
       
def transformelg(a,b,Ordre):
    x,w=lglnodes(Ordre)
    x_transformed = (b - a) / 2 * x + (a + b) / 2
    w_transformed = (b - a) / 2 * w
    return x_transformed, w_transformed



mat=np.eye(1)
rho=np.eye(1)
c=1
dt=0.01
dx=1
Longueur=10
L=int(Longueur/dx) #nombre d'élements mais L+1 noeuds
r=3
Duree=30
N=int(Duree/dt)

f=0.5
periode=1/f
x=Longueur/2
exemplepoint=5
amplitude=100
diracIndice=int((x/dx)+(r-1)*((x/dx)-1))



#discrétisation d'espace respectant les noeuds
Mesh=[0]
X=[] #liste des coordonées de chaque point par éléments
N_point=L+1+(r-1)*(L)
Id_noeud=[i*dx for i in range(0,N_point,r)]


for i in range(0,L):
    x,w=transformelg(i*dx,(i+1)*dx,r)
    X.append(x)
    for t in range(1,r+1):
        Mesh.append(x[t])
    

# pour=[]



# K_local=[] # stocke les K pour chaque éléments 
# M_local=[] #stock les M pour chaque élements
# for t in range(len(X)):
#     r=len(X[0])
#     ktemp=np.zeros((r,r))
#     mtemps=np.zeros((r,r))
#     coords=np.array(X[t]).reshape(r,1)
#     for i in range(r):
#         c=np.array(dphi[i,:]).reshape(r,1) #vecteur dphi
#         d=np.array(phi[:,i]).reshape(r,1) #vecteur phi
#         J=c.T@coords             #Jacobien
#         aux = np.linalg.inv(J)@c.T
#         baux=np.linalg.inv(J)@d.T
#         ktemp=ktemp+(w[i]*np.linalg.det(J)*aux.T@mat@aux)
#         mtemps=mtemps+(w[i]*np.linalg.det(J)* baux.T@rho@baux)
#     K_local.append(ktemp)
#     M_local.append(mtemps)
#     pour.append(J)

# r=len(K_local[0][:,0])-1
# mat=np.zeros((N_point,N_point))
# for i in range(len(K_local)):
#     mat[r*i:r+1+r*i,r*i:r+1+r*i]=K_local[i][:,:]
#     if i<len(K_local)-1:
#          mat[r+r*i,r+r*i]=K_local[i][r,r]+K_local[i+1][0,0]
# mata=np.zeros((N_point,N_point))
# for i in range(len(M_local)):
#     mata[r*i:r+1+r*i,r*i:r+1+r*i]=M_local[i][:,:]
#     if i<len(M_local)-1:
#          mata[r+r*i,r+r*i]=M_local[i][r,r]+M_local[i+1][0,0]





def fabricmatrice(X,phi,dphi,rho,mat):
    K_local=[] # stocke les K pour chaque éléments 
    M_local=[] #stock les M pour chaque élements
    r=len(X[0])
    for t in range(len(X)):
        ktemp=np.zeros((r,r))
        mtemps=np.zeros((r,r))
        coords=np.array(X[t]).reshape(r,1)
        for i in range(r):
            c=np.array(dphi[i,:]).reshape(r,1)
            d=np.array(phi[:,i]).reshape(r,1)
            J=c.T@coords #Jacobienne
            aux = np.linalg.inv(J)@ c.T
            baux=np.linalg.inv(J)@d.T
            ktemp=ktemp+(w[i]*np.linalg.det(J)*aux.T@mat@aux)
            mtemps=mtemps +(w[i]*np.linalg.det(J)* baux.T@rho@baux)
        K_local.append(ktemp)
        M_local.append(mtemps)
    return K_local,M_local

   
def assemble(Y,taille):
    r=len(Y[0][:,0])-1
    mat=np.zeros((taille,taille))
    for i in range(len(Y)):
        mat[r*i:r+1+r*i,r*i:r+1+r*i]=Y[i][:,:]
    for i in range(len(Y)):
        if i<len(Y)-1:
            mat[r+r*i,r+r*i]=Y[i][r,r]+Y[i+1][0,0]
    return mat



K_local,M_local=fabricmatrice(X,phi,dphi,rho,mat)   
K=assemble(K_local,N_point) #Matrice de rigidité
M=assemble(M_local,N_point)
M_inverse=np.diag(1/np.diag(M))

U=np.zeros((N_point,N))


#etablissement de la matrice contenant le sinus porte
dirac=np.array([0 for i in range(N_point)])
dirac[diracIndice]=1
dirac=dirac.reshape(N_point,1)
t = np.linspace(0, Duree, N)
F=amplitude*np.sin(2* np.pi * t / periode) * (t <= periode)
F=F.reshape(1,N)
DiracF=np.dot(dirac,F)

for i in range(2,N):
    U[:,i]=dt**2*np.dot(M_inverse,DiracF[:,i-1]-c*K@U[:,i-1]) +2*U[:,i-1]-U[:,i-2]
    
# Observation en un point quelconque 

plt.figure()

plt.plot(t,U[int(exemplepoint/dx),:],"b",label="amplitude en x=")
plt.xlabel("temps")
plt.ylabel("déplacement")

plt.legend()
plt.show()
plt.figure()
plt.plot(t,DiracF[diracIndice,:],"r",label="sinus porte")
plt.xlabel("temps")
plt.ylabel("dépmacement")
plt.title("probleme 1D propagation onde")
plt.legend()
plt.show()    