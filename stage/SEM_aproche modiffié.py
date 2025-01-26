# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 10:56:08 2024

@author: Utilisateur
"""
import math
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import solution_exacte as exact




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

def Force(x,L,duree,f,amplitude,r,dt,dx,t):
    N_element=int(L/dx)
    lvecteur=N_element*r+1
    cvecteur=int(duree/dt)
    periode=1/f
    diracIndice=int(x/dx)*r+1
    # etablissement de la matrice contenant le sinus porte
    dirac=np.array([0 for i in range(lvecteur)])
    dirac[diracIndice]=1
    dirac=dirac.reshape(lvecteur,1)
    
    F=amplitude*np.sin(2* np.pi * t / periode) * (t <= periode)
    F=F.reshape(1,cvecteur+1)
    DiracF=np.dot(dirac,F)
    
    return DiracF
def maillage(L,r,dx):
    Mesh=[0] 
    X=[]    #liste des coordonées de chaque point par éléments
    for i in range(0,L):
        x,w=transformelg(i*dx,(i+1)*dx,r)
        X.append(x)
        for t in range(1,r+1):
            Mesh.append(x[t])
    return Mesh,X,w

def approcheModifié(U,K,M_inverse,F,dt,c,k):
    c_2=c**2
    Dernier=len(U[:,0])-1
    N=M_inverse@K
    for i in range(2,len(U[0,:])):
        U[:,i]=(dt**2)*((M_inverse@F[:,i-1])-(c_2*N@U[:,i-1])) +2*U[:,i-1]-U[:,i-2]
        if k>1: 
            for j in range(2,k+1):
                coef=2*((-1)*c_2)**j*dt**(2*j)/math.factorial(2*j)
                U[:,i]+=coef*np.linalg.matrix_power(N,j)@U[:,i-1]
        U[0,i]=0
        U[Dernier,i]=0
    return U

def Matrice_deplacement(c,rho,mat,duree,Longueur,r,dt,dx,f,a,amplitude,t,k):
    N=int(duree/dt)
    L=int(Longueur/dx)
    N_point=L*r+1
    x,w= lglnodes(r);
    phi,dphi =ShapeFunctionLagrange(x,method_type='SEM')
    Mesh,X,w=maillage(L,r,dx)
    K_local,M_local=fabricmatrice(X,w,phi,dphi,rho,mat)   
    K=assemble(K_local,N_point) #Matrice de rigidité
    M=assemble(M_local,N_point)
    M_inverse=np.diag(1/np.diag(M))
    U=np.zeros((N_point,N+1))
    F=Force(a,Longueur,duree,f,amplitude,r,dt,dx,t)
    U=approcheModifié(U,K,M_inverse,F,dt,c,k)
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
amplitude=100/22
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

# F=Force(a,Longueur,Duree,f,amplitude,r,dt,dx,t)
# # U,Mesh=Matrice_deplacement(c,rho,mat,Duree,Longueur,r,dt,dx,f,a,amplitude,t,2)

# U=approcheModifié(U, K, M_inverse, F, dt, c, 2)


def erreur_MEA(K,M_inverse,c,duree,Longueur,rho,mat,f,A,a,dt,dx,r,t,k):
    omega=2*np.pi*f
    T=1/f
    N=int(duree/dt)
    N_point=int(Longueur/dx)*r+1
    U=np.zeros((N_point,N+1))
    F=Force(a,Longueur,duree,f,A,r,dt,dx,t)
    U=approcheModifié(U,K,M_inverse,F,dt,c,k)
    U_exact=exact.Matrice_deplacement(Mesh,a,c,L,A,omega,T,t)
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
    plt.figure()
    plt.scatter(alpha,Y[0],label="ordre=2",marker='x')
    for i in range(1,len(Y)):
        plt.scatter(alpha,Y[i],label=f'ordre={i*2}',marker='x')
    plt.xlabel("alpha")
    plt.ylabel("Log(norme_L2(erreur)")
    # plt.ylim(0,40)
    plt.grid()
    plt.legend()
    plt.show()
    plt.figure()
    plt.scatter(temps,Y[0],label="ordre=2",marker='x')
    for i in range(1,len(Y)):
        plt.scatter(temps,Y[i],label=f'ordre={i*2}',marker='x')
    plt.xlabel("dt")
    plt.ylabel("Log(norme_L2(erreur)")
    # plt.ylim(0,40)
    plt.grid()
    plt.legend()
    plt.show()
    
    return
    
    
def select(L,seuil):
    i=0
    while L[i]<seuil and i<len(L)-1:
        i=i+1
    for j in range(i,len(L)-1):
        del L[i]
    return i

temps,alpha,Y=test_MEA (c,Duree,Longueur,rho,mat,f,amplitude,a,0.05,2,5,100)

affiche_1(alpha,temps,Y)


# seuil=50

 


# plt.figure()
# plt.scatter([1,2,3,4,5], [k0,k2,k3,k4,k5], color='red')
# plt.annotate("k0", (1, k0), textcoords="offset points", xytext=(0, 10), ha='center')
# plt.annotate("k2", (2, k2), textcoords="offset points", xytext=(0, 10), ha='center')
# plt.annotate("k3", (3, k3), textcoords="offset points", xytext=(0, 10), ha='center')
# plt.annotate("k4", (4, k4), textcoords="offset points", xytext=(0, 10), ha='center')
# plt.annotate("k5", (5, k5), textcoords="offset points", xytext=(0, 10), ha='center')
# plt.show()


# # Visualisation dU déplacement sur tout le
# plt.figure(figsize=(10, 6))
# plt.title(f"Solution élements spectraux de l'équation de d'Alembert r={r} ")
# plt.pcolormesh(t, Mesh, U, shading='auto', cmap='jet')
# plt.colorbar(label="Amplitude")
# plt.xlabel("Temps (s)")
# plt.ylabel("Position (x)")
# plt.grid(True)
# plt.show()





# #Affichage de la source, déplacement en 1 point


# plt.figure()
# plt.plot(t,U[int(exemplepoint/dx),:],label=f"amplitude en x={exemplepoint}")
# plt.xlabel("temps")
# plt.ylabel("déplacement")
# plt.grid()
# plt.legend()
# plt.show()
# plt.figure()
# plt.plot(t,DiracF[diracIndice,:],"r--",label="sinus porte")
# plt.xlabel("temps")
# plt.ylabel("dépmacement")
# plt.title("probleme 1D propagation onde")
# plt.legend()
# plt.show() 

# fig, ax=plt.subplots(subplot_kw={"projection":"3d"})
# x,y=np.meshgrid(t,Mesh)
# surf=ax.plot_surface(x,y,U)
# plt.show()


# fig=plt.figure()
# plt.xlabel("longueur")
# plt.ylabel("deplacement")
# plt.title("probleme 1D propagation onde")
# line,=plt.plot([],[])
# plt.xlim(0,Longueur)
# plt.ylim(-amplitude*2,amplitude*2)

# def animate(i):
#      line.set_data(Mesh,U[:,i])
#      return line,
# ani = animation.FuncAnimation(fig, animate, frames=range(0, N, 100), interval=10, blit=True, repeat=False)

# # Enregistrement de l'animation
# ani.save("test_dt=0.1.mp4", writer="ffmpeg", fps=30)   