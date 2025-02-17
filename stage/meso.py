# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 13:19:29 2024

@author: Utilisateur
"""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pickle 
import solution_exact_bis as exact




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

def transformelg(a,b,Ordre):
    x,w=lglnodes(Ordre)
    x_transformed = (b - a) / 2 * x + (a + b) / 2
    w_transformed = (b - a) / 2 * w
    return x_transformed, w_transformed

def fabricmatrice(X,w,phi,dphi,rho,mat):
    K_local=[] # stocke les K pour chaque éléments 
    M_local=[] #stock les M pour chaque élements
    r=len(X[0])
    for t in range(len(X)):
        ktemp=np.zeros((r,r))
        mtemps=np.zeros((r,r))
        coords=np.array(X[t]).reshape(r,1)
        for i in range(r):
            c=np.array(dphi[i,:]).reshape(r,1)
            d=np.array(phi[i,:]).reshape(r,1)
            J=c.T@coords #Jacobien
            aux = np.linalg.inv(J)@ c.T
            ktemp=ktemp+(w[i]*np.linalg.det(J)*aux.T@mat@aux)
            mtemps += w[i] * np.linalg.det(J) * d @ d.T
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


    

def maillage(L,r,dx):
    Mesh=[0] 
    X=[]    #liste des coordonées de chaque point par éléments
    for i in range(0,L):
        x,w=transformelg(i*dx,(i+1)*dx,r)
        X.append(x)
        for t in range(1,r+1):
            Mesh.append(x[t])
    return Mesh,X,w

def approcheModifié(duree,Longueur,r,f,a,ua,K,M_inverse,dt,dx,c,k):
    c_2=c**2
    N_temps=int(duree/dt)
    L=int(Longueur/dx)
    N_espace=L*r+1
    a_=int(a/dx)*r+1
    U=np.zeros((N_espace,N_temps+1))
    N=M_inverse@K
    if k>1:
        liste=[]
        for j in range(2,k+1):
            coef=2*((-1)*c_2)**j*dt**(2*j)/math.factorial(2*j)
            liste.append(coef*np.linalg.matrix_power(N,j))
        for i in range(2,N_temps):
            U[:,i]=-(dt**2)*(c_2)*N@U[:,i-1] +2*U[:,i-1]-U[:,i-2]
            for j in range(k-2):
                U[:,i] += liste[j]@U[:,i-1]
            U[a_,i]=ua[i]
    else:
        for i in range(2,N_temps):
            U[:,i]=-(dt**2)*(c_2)*N@U[:,i-1] +2*U[:,i-1]-U[:,i-2]
            U[a_,i]=ua[i]
        
    return U

def condition(t,A,f):
    omega=f*2*np.pi
    T=1/f
    valeur=A*np.sin(omega*t) *(t<=T)
    return valeur

def Matrice_deplacement(c,rho,mat,duree,Longueur,r,dt,dx,f,a,amplitude,t,k):
    N=int(duree/dt)
    L=int(Longueur/dx)
    N_point=L*r+1
    ua=[condition(i*dt,amplitude,f) for i in range(N+1)]
    x,w= lglnodes(r);
    phi,dphi =ShapeFunctionLagrange(x,method_type='SEM')
    Mesh,X,w=maillage(L,r,dx)
    K_local,M_local=fabricmatrice(X,w,phi,dphi,rho,mat)   
    K=assemble(K_local,N_point) #Matrice de rigidité
    M=assemble(M_local,N_point)
    M_inverse=np.diag(1/np.diag(M))
    U=np.zeros((N_point,N+1))
    U=approcheModifié(duree,Longueur,r,f,a,ua,K,M_inverse,dt,dx,c,k)
    return U,Mesh
               






#Paramètre simul/matériaux/Onde/Source/
dt=0.01
dx=0.05
r=2

mat=np.eye(1)
rho=np.eye(1)
Longueur=10
L=int(Longueur/dx) #nombre d'élements mais L+1 noeuds
Duree=30
N=int(Duree/dt) #nombre d'intervale de temps
c=1
N_point=L*r+1  #nombre de point

f=0.5
periode=1/f
a=Longueur/2
exemplepoint=1
amplitude=100
diracIndice=int((a/dx)*r+1)

#discrétisation du temps
t = np.linspace(0, Duree, N+1)
#discrétisation en espace
Mesh=maillage(L,r,dx)[0]


x,w= lglnodes(r);
phi,dphi =ShapeFunctionLagrange(x,method_type='SEM')
Mesh,X,w=maillage(L,r,dx)
K_local,M_local=fabricmatrice(X,w,phi,dphi,rho,mat)   
K=assemble(K_local,N_point) #Matrice de rigidité
M=assemble(M_local,N_point)
M_inverse=np.diag(1/np.diag(M))

U=np.zeros((N_point,N+1))


# # U,Mesh=Matrice_deplacement(c,rho,mat,Duree,Longueur,r,dt,dx,f,a,amplitude,t,2)




def erreur_MEA(K,M_inverse,c,duree,Longueur,rho,mat,f,A,a,dt,dx,r,t,k):
    omega=2*np.pi*f
    T=1/f
    N=int(duree/dt)
    N_point=int(Longueur/dx)*r+1
    ua=[condition(i*dt,A,f) for i in range(N)]
    U=np.zeros((N_point,N+1))
    U=approcheModifié(duree,Longueur,r,f,a,ua,K,M_inverse,dt,dx,c,k)
    U_exact=exact.matrice_depl(t,Mesh,a,c,A,omega,T,L)
    return math.log(np.linalg.norm(U_exact-U) ,10)
    
# print(erreur_MEA(K,M_inverse,c,Duree,Longueur,rho,mat,f,amplitude,a,dt,dx,r,t,2))

    

        
    




def test_MEA (c,Duree,Longueur,rho,mat,f,amplitude,a,dx,r,k_max,nombre_point):
    temps=[0.01+i*0.001 for i in range(0,nombre_point)]
   
    alpha=[c*dt/dx for dt in temps]
    Y=[]
    for i in range(1,k_max+1):
        MEA=[]
        for dt in temps: 
            Nb=int(Duree/dt)
            t=np.linspace(0,Duree,Nb+1)
            MEA.append(erreur_MEA(K,M_inverse,c,Duree,Longueur,rho,mat,f,amplitude,a,dt,dx,r,t,i) )
        
        Y.append(MEA)
        MEA=[]
    return temps,alpha,Y

def affiche_1(alpha,temps,Y):
    fig, ax = plt.subplots()
    plt.scatter(alpha,Y[0],label="ordre=2",marker='x')
    for i in range(1,len(Y)):
        plt.scatter(alpha,Y[i],label=f'ordre={i*2+2}',marker='x')
    plt.xlabel("alpha")
    plt.ylabel("Log(norme_L2(erreur)")
    plt.grid()
    plt.legend()
    plt.title("Erreur en fonction du CFL")
    plt.savefig("Stabilité_alpha",dpi=300)
    with open('Stabilité_alpha.pkl', 'wb') as f:
        pickle.dump(fig, f)
    plt.close()
    fig, ax = plt.subplots()
    plt.scatter(temps,Y[0],label="ordre=2",marker='x')
    for i in range(1,len(Y)):
        plt.scatter(temps,Y[i],label=f'ordre={i*2+2}',marker='x')
    plt.xlabel("dt (en s) ")
    plt.ylabel("Log(norme_L2(erreur)")
    plt.grid()
    plt.legend()
    plt.title("Erreur en fonction du pas de temps")
    plt.savefig("Stabilité_dt",dpi=300)
    plt.savefig("Stabilité_alpha",dpi=300)
    with open('Stabilité_dt.pkl', 'wb') as f:
        pickle.dump(fig, f)
    plt.close()
    
    return

    
def select(L,seuil):
    i=0
    while L[i]<seuil and i<len(L)-1:
        i=i+1
    for j in range(i,len(L)-1):
        del L[i]
    return i

temps,alpha,Y=test_MEA (c,Duree,Longueur,rho,mat,f,amplitude,a,dx,r,10,300)

affiche_1(alpha,temps,Y)


# seuil=50