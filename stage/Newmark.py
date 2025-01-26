import numpy as np
import matplotlib.pyplot as plt


c=1
dt=0.001
dx=0.1
alpha=c*dt/dx
print("alpha=",alpha)
Longueur=30
duree=30

f=0.5
periode=1/f
x=Longueur/2
exemplepoint=1
amplitude=100
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
acc_prec=np.zeros((L,1))
v=np.zeros((L,1))

# resolution avec algo de newmark gamma=1/2 beta=0
for n in range(1,N):
    acc_prec[:,0]=acc[:,0]
    acc[:,0]=DiracF[:,n-1]-np.dot(K,U[:,n-1])
    v[:,0],U[:,n]=v[:,0]+dt*0.5*(acc[:,0]+acc_prec[:,0]), U[:,n-1]+dt*v[:,0]+(dt**2)*0.5*acc_prec[:,0]
    

X=[i*dx for i in range(L)]
plt.figure()

plt.plot(t,U[int(exemplepoint//dx),:],"b",label="amplitude en x")
plt.xlabel("temps")
plt.ylabel("déplacement")
plt.title("probleme 1D propagation onde Newmark")
plt.legend()
plt.show()
plt.figure()
plt.plot(t,DiracF[diracIndice,:],"r",label="sinus porte")
plt.xlabel("temps")
plt.ylabel("force")
plt.title("probleme 1D propagation onde Newmark")
plt.legend()
plt.show() 

# fig=plt.figure()
# plt.xlabel("longueur")
# plt.ylabel("deplacement")
# plt.title("probleme 1D propagation onde")
# line,=plt.plot([],[])
# plt.xlim(0,Longueur)
# plt.ylim(-amplitude*0.01,amplitude*0.01)


# def animate(i):
#      line.set_data(X,U[:,i])
#      return line,
# ani = animation.FuncAnimation(fig, animate, frames=range(0, N, 100), interval=10, blit=True, repeat=False)

# # Enregistrement de l'animation
# ani.save("propagation_ondes_Newmark.mp4", writer="ffmpeg", fps=30)
# plt.show()


fig, ax=plt.subplots(subplot_kw={"projection":"3d"})
x,y=np.meshgrid(t,X)
surf=ax.plot_surface(x,y,U)
plt.show()   
