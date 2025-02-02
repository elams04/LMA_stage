# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 08:41:46 2024

@author: Utilisateur
"""


import numpy as np



import matplotlib

# Réinitialiser les paramètres à leurs valeurs par défaut
matplotlib.rc_file_defaults()
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




       
def transformelg(a,b,Ordre):
    x,w=lglnodes(Ordre)
    x_transformed = (b - a) / 2 * x + (a + b) / 2
    w_transformed = (b - a) / 2 * w
    return x_transformed, w_transformed

def predicteur_corecteur(M_inverse,K,F,c,dt,U):
    L=len(U[:,0])
    N=len(U[0,:])
    acc=np.zeros((L,1))
    v_pred=np.zeros((L,1))
    v=np.zeros((L,1))
    u_pred=np.zeros((L,1))
    gamma=1/2
    beta=0

    #resolution avec algo de newmark
    for n in range(2,N):    
        v_pred[:,0], u_pred[:,0]=v[:,0]+dt*(1-gamma)*acc[:,0], U[:,n-1]+dt*v[:,0]+(dt**2)*0.5*(1-2*beta)*acc[:,0]
        
        #calcul accéleration
        acc[:,0]=M_inverse@(F[:,n-1]-np.dot(c*K,u_pred[:,0]))
        
        #corections
        U[:,n]=(beta*dt**2)*acc[:,0]+u_pred[:,0]
        v[:,0]=v_pred[:,0]+gamma*dt*acc[:,0]
    return U

def condition(t,A,f):
    omega=f*2*np.pi
    T=1/f
    valeur=A*np.sin(omega*t) *(t<=T)
    return valeur

def FDM(duree,Longueur,r,f,a,ua,K,M_inverse,dt,dx,c):
    N=int(duree/dt)
    L=int(Longueur/dx)
    N_point=L*r+1
    a_=int(a/dx)*r
    U=np.zeros((N_point,N+1))
    for i in range(2,N+1):
        if i*dt<=(1/f) :
            U[:,i]=-(dt**2)*(c**2)*M_inverse@K@U[:,i-1] +2*U[:,i-1]-U[:,i-2]
            U[a_,i]=ua[i]
        else:
            U[:,i]=-(dt**2)*(c**2)*M_inverse@K@U[:,i-1] +2*U[:,i-1]-U[:,i-2]
        
        
        
    return U


def maillage(L,r,dx):
    Mesh=[0] 
    X=[]    #liste des coordonées de chaque point par éléments
    for i in range(0,L):
        x,w=transformelg(i*dx,(i+1)*dx,r)
        X.append(x)
        for t in range(1,r+1):
            Mesh.append(x[t])
    return Mesh,X,w



def fabricmatrice(X,w,phi,dphi,rho,mat):
    K_local=[] # stocke les K pour chaque éléments 
    M_local=[] #stock les M pour chaque élements
    r=len(X[0])
    for t in range(len(X)):
        ktemp=np.zeros((r,r))
        mtemp=np.zeros((r,r))
        coords=np.array(X[t]).reshape(r,1)
        for i in range(r):
            c=np.array(dphi[i,:]).reshape(r,1)
            d=np.array(phi[i,:]).reshape(r,1)
            J=c.T@coords #Jacobien
            aux = np.linalg.inv(J)@ c.T
            ktemp=ktemp+(w[i]*np.linalg.det(J)*aux.T@mat@aux)
            mtemp += w[i] * np.linalg.det(J) * d @rho@ d.T
        K_local.append(ktemp)
        M_local.append(mtemp)
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


def Matrice_deplacement(c,rho,mat,duree,Longueur,r,dt,dx,f,a,amplitude,t):
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
    U=FDM(Duree,Longueur,r,a,ua,K,M_inverse,dt,dx,c)

    return U,Mesh

if __name__=="__main__":
    #Paramètre simul/matériaux/Onde/Source/
    nombre_point_minimun=200
    
    
    #Soit on fixe c et don on laisse mat et rho =1 ou on fixe rho et mat et donc c=1
    mat=np.eye(1)
    rho=np.eye(1)
    c=1
    
    Longueur=10
    dt=0.01
    dx=0.1
    r=2
    
    L=int(Longueur/dx) #nombre d'élements mais L+1 noeuds
    Duree=30
    N=int(Duree/dt) #nombre d'intervale de temps
    
    
    alpha=c*dt/dx
    

    f=0.5
    periode=1/f
    a=Longueur/2    #endroit où de la source
    exemplepoint=1  # observation en un point quelconque
    instant=20       #observation en un temps quelconque
    amplitude=10


    #discrétisation du temps
    t = np.linspace(0, Duree, N+1)

    # #discrétisation d'espace respectant les noeuds
    # Mesh=[0]
    # X=[]
    N_point=L*r+1
    if N_point<nombre_point_minimun:
        print("Attention Pas assez de point")
    # Id_noeud=[i*dx for i in range(0,N_point,r)]


                   

                   

    x,w= lglnodes(r);
    phi,dphi =ShapeFunctionLagrange(x,method_type='SEM')
    Mesh,X,w=maillage(L,r,dx)
    K_local,M_local=fabricmatrice(X,w,phi,dphi,rho,mat)   
    K=assemble(K_local,N_point) #Matrice de rigidité
    M=assemble(M_local,N_point)
    M_inverse=np.diag(1/np.diag(M))
    ua=[condition(i*dt,amplitude,f) for i in range(N+1)]

    U=FDM(Duree,Longueur,r,f,a,ua,K,M_inverse,dt,dx,c)





    # #Algo predicteur correcteur 
    # U=predicteur_corecteur(M_inverse,K,DiracF,c,dt,U)

    # Algo difference finie


    # U,Mesh=Matrice_deplacement(c,rho,mat,Duree,Longueur,r,dt,dx,f,a,amplitude,t)


    # Observation en un point quelconque 

    
    plt.figure(figsize=(8, 5))
    plt.plot(t,U[int(exemplepoint/dx)*r+1,:],"b",label=f'amplitude en x={exemplepoint}')
    plt.xlabel("temps")
    plt.ylabel("déplacement")
    plt.title(f'Solution SEM en x={exemplepoint}')
    plt.grid(True)
    plt.legend()
    plt.show()
    
    # Observation en un temps quelconque 

    
    plt.figure(figsize=(8, 5))
    plt.plot(Mesh,U[:,int(instant/dt)],"b",label=f'corde en t={instant}')
    plt.xlabel("temps")
    plt.ylabel("déplacement")
    plt.title(f'Solution SEM en x={exemplepoint}')
    plt.grid(True)
    plt.legend()
    plt.show()
    # plt.figure()
    # plt.plot(t,DiracF[diracIndice,:],"r",label="sinus porte")
    # plt.xlabel("temps")
    # plt.ylabel("dépmacement")
    # plt.title("probleme 1D propagation onde")
    # plt.legend()
    # plt.show()
     
    # Visualisation dU déplacement sur tout le
    plt.figure(figsize=(10, 6))
    plt.title(f"Solution élements spectraux de l'équation de d'Alembert r={r} ")
    plt.pcolormesh(t, Mesh, U, shading='auto', cmap='jet')
    plt.colorbar(label="Amplitude")
    plt.xlabel("Temps (s)")
    plt.ylabel("Position (x)")
    plt.grid(True)
    plt.show()

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
    # plt.ylim(-amplitude*5,amplitude*5)


    # def animate(i):
    #      line.set_data(Mesh,U[:,i])
    #      return line,
    # ani = animation.FuncAnimation(fig, animate, frames=range(0, N, 10), interval=10, blit=True, repeat=False)

    # # Enregistrement de l'animation
    # ani.save("propagation_ondes_SEM_F.mp4", writer="ffmpeg", fps=30)   