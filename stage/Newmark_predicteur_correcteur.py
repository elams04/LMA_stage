# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 08:58:42 2024

@author: Utilisateur
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

c=1
dt=0.0010
dx=0.01
alpha=c*dt/dx
print("alpha=",alpha)
Longueur=10
duree=30

f=0.5
periode=1/f
x=Longueur/2
exemplepoint=1
amplitude=10
diracIndice=int(x//dx)
L=int(Longueur/dx)
N=int(duree/dt)

# etablissement de la matrice K issue de la discrétisation d'ordre 2 d'espace
diagoinf=np.ones(L-1)
diagoprincipale=-2*np.ones(L)
K=-(c/dx)**2*(np.diag(diagoprincipale)+np.diag(diagoinf,-1)+np.diag(diagoinf,1))
K[0,0]=0
K[0,1]=0
K[L-1,L-1]=0
K[L-1,L-2]=0

# etablissement de la matrice contenant la force en sinus porte
dirac=np.array([0 for i in range(L)])
dirac[diracIndice]=1
dirac=dirac.reshape(L,1)
t = np.linspace(0, duree, N)
F=amplitude*np.sin(2* np.pi * t / periode) * (t <= periode)
F=F.reshape(1,N)
DiracF=np.dot(dirac,F)

U=np.zeros((L,N))
acc=np.zeros((L,1))
v_pred=np.zeros((L,1))
v=np.zeros((L,1))
u_pred=np.zeros((L,1))
gamma=1/2
beta=0
def condition(t,A,f):
    omega=f*2*np.pi
    T=1/f
    valeur=A*np.sin(omega*t) *(t<=T)
    return valeur


ua=[condition(i*dt,amplitude,f) for i in range(int(duree/dt)+1)]
#resolution avec algo de newmark
for n in range(1,N):    
    v_pred[:,0], u_pred[:,0]=v[:,0]+dt*(1-gamma)*acc[:,0], U[:,n-1]+dt*v[:,0]+(dt**2)*0.5*(1-2*beta)*acc[:,0]
    
    #calcul accéleration
    acc[:,0]=-np.dot(K,u_pred[:,0])
    
    #corections
    U[:,n]=(beta*dt**2)*acc[:,0]+u_pred[:,0]
    v[:,0]=v_pred[:,0]+gamma*dt*acc[:,0]
    U[diracIndice,n]=ua[n]


X=[i*dx for i in range(L)]
plt.figure()

plt.plot(t,U[int(exemplepoint//dx),:],"b",label="amplitude en x")
plt.xlabel("temps")
plt.ylabel("déplacement")
plt.title("probleme 1D propagation onde Newmark prédicteur correcteur")
plt.grid(True)
plt.legend()
plt.show()
# plt.figure()
# plt.plot(t,DiracF[diracIndice,:],"r",label="sinus porte")
# plt.xlabel("temps")
# plt.ylabel("force")
# plt.title("probleme 1D propagation onde Newmark")
# plt.legend()
# plt.show() 

# fig, ax=plt.subplots(subplot_kw={"projection":"3d"})
# x,y=np.meshgrid(t,X)
# surf=ax.plot_surface(x,y,U)
# plt.show()

# fig=plt.figure()
# plt.xlabel("longueur")
# plt.ylabel("deplacement")
# plt.title("probleme 1D propagation onde")
# line,=plt.plot([],[])
# plt.xlim(0,Longueur)
# plt.ylim(-amplitude*0.1,amplitude*0.1)
# def animate(i):
#      line.set_data(X,U[:,i])
#      return line,
# ani = animation.FuncAnimation(fig, animate, frames=range(0, N, 100), interval=10, blit=True, repeat=False)

# # Enregistrement de l'animation
# ani.save("propagation_ondes_Newmark_predicteur_corecteur.mp4", writer="ffmpeg", fps=30)
# plt.show()