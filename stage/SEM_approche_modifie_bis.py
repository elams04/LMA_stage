# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 13:01:44 2024

@author: Utilisateur
"""

import math
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
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
            mtemps += w[i] * np.linalg.det(J) * d @rho@ d.T
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
    a_=int(a/dx)*r
    U=np.zeros((N_espace,N_temps+1))
    N=M_inverse@K
    if k>1:
        liste=[]
        for j in range(2,k+1):
            coef=2*(((-1)*c_2)**j)*(dt**(2*j))/math.factorial(2*j)
            liste.append(coef*np.linalg.matrix_power(N,j))
        for i in range(2,N_temps+1):
            if i*dt<=(1/f):
                U[:,i]=-(dt**2)*c_2*N@U[:,i-1] +2*U[:,i-1]-U[:,i-2]
                for j in range(len(liste)):
                    U[:,i] += liste[j]@U[:,i-1]
                U[a_,i]=ua[i]
            else:
                U[:,i]=-(dt**2)*(c_2)*N@U[:,i-1] +2*U[:,i-1]-U[:,i-2]
                for j in range(len(liste)):
                    U[:,i] += liste[j]@U[:,i-1]
            
            
    else:
        for i in range(2,N_temps+1):
            if i*dt<=(1/f):
                U[:,i]=-(dt**2)*(c_2)*N@U[:,i-1] +2*U[:,i-1]-U[:,i-2]
                U[a_,i]=ua[i]
            else:
                U[:,i]=-(dt**2)*(c_2)*N@U[:,i-1] +2*U[:,i-1]-U[:,i-2]
            
        
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
               






if __name__=="__main__":
    #Paramètre simul/matériaux/Onde/Source/
    dt=0.001
    dx=0.02
    r=2
    k=2

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
    instant=20
    amplitude=10

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
    ua=[condition(i*dt,amplitude,f) for i in range(N+1)]
    U=approcheModifié(Duree, Longueur, r, f, a, ua, K, M_inverse, dt, dx, c, k)

    # U,Mesh=Matrice_deplacement(c,rho,mat,Duree,Longueur,r,dt,dx,f,a,amplitude,t,2)




    # Visualisation dU déplacement sur tout le
    plt.figure(figsize=(10, 6))
    plt.title(f"Solution élements spectraux de l'équation de d'Alembert r={r} ")
    plt.pcolormesh(t, Mesh, U, shading='auto', cmap='jet')
    plt.colorbar(label="Amplitude")
    plt.xlabel("Temps (s)")
    plt.ylabel("Position (x)")
    plt.grid(True)
    plt.show()





    # #Affichage  déplacement en 1 point


    plt.figure(figsize=(8, 5))
    plt.plot(t,U[int(exemplepoint/dx)*r,:],label=f"amplitude en x={exemplepoint}")
    plt.xlabel("temps")
    plt.ylabel("déplacement")
    plt.title(f'Solution SEM r={r} ordre={k}')
    plt.grid()
    plt.legend()
    plt.show()
    # Affichage en un temps quelconque
    plt.figure(figsize=(8, 5))
    plt.plot(Mesh,U[:,int(instant/dt)],"b",label=f'corde en t={instant}')
    plt.xlabel("temps")
    plt.ylabel("déplacement")
    plt.title(f'Solution SEM k={k} en x={exemplepoint}')
    plt.grid(True)
    plt.legend()
    plt.show()


    fig, ax=plt.subplots(subplot_kw={"projection":"3d"})
    x,y=np.meshgrid(t,Mesh)
    surf=ax.plot_surface(x,y,U)
    plt.show()


    fig=plt.figure()
    plt.xlabel("longueur")
    plt.ylabel("deplacement")
    plt.title("probleme 1D propagation onde")
    line,=plt.plot([],[])
    plt.xlim(0,Longueur)
    plt.ylim(-amplitude*2,amplitude*2)

    def animate(i):
         line.set_data(Mesh,U[:,i])
         return line,
    ani = animation.FuncAnimation(fig, animate, frames=range(0, N, 50), interval=0.01, blit=True, repeat=False)

    # Enregistrement de l'animation
    ani.save("propagation_ondesapproche modifé k={k}.mp4", writer="ffmpeg", fps=30)
    plt.show()