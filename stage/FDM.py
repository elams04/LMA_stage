import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


c=1
pastemps=0.01
pasespace=0.1
alpha=c*pastemps/pasespace
print("alpha=",alpha)
Longueur=10
duree=30

f=0.5
periode=1/f
x=Longueur/2
exemplepoint=1
amplitude=100
diracIndice=int(x/pasespace)
lvecteur=int(Longueur/pasespace)
cvecteur=int(duree/pastemps)


# # Conditions de bords
# a=0
# b=0

# #etablissement de la matrice A
# gamma=alpha**2
# diagoinf=gamma*np.ones(lvecteur)
# beta=beta = 2*(1-alpha**2)
# diagoprincipale=beta*np.ones(lvecteur+1)
# A=(np.diag(diagoprincipale)+np.diag(diagoinf,-1)+np.diag(diagoinf,1))
# A[0,1]=0
# A[lvecteur,lvecteur-1]=0
# A[0,0]=a
# A[lvecteur,lvecteur]=b




# # etablissement de la matrice contenant le sinus porte
# dirac=np.array([0 for i in range(lvecteur+1)])
# dirac[diracIndice]=1/pasespace
# dirac=dirac.reshape(lvecteur+1,1)
# t = np.linspace(0, duree, cvecteur+1)
# F=amplitude*np.sin(2* np.pi * t / periode) * (t <= periode)
# F=F.reshape(1,cvecteur+1)
# DiracF=np.dot(dirac,F)


# u0=0 pour tout x v0=0 pour tout x à t=0

# U=np.zeros((lvecteur+1,cvecteur+1))

# X=[i*pasespace for i in range(lvecteur+1)]

# #Resolution du prblm

# for i in range(2,cvecteur+1):
#     U[:,i] = A @ U[:,i-1] - U[:,i-2] + (pastemps**2) * DiracF[:,i-1]
#     U[0,i]=0
#     U[lvecteur,i]=0

def Matrice_deplacement(c,duree, L,DiracF,dt,dx):
    alpha=c*dt/dx
    lvecteur=int(L/dx)
    cvecteur=int(duree/dt)
    #etablissement de la matrice A
    gamma=alpha**2
    diagoinf=gamma*np.ones(lvecteur)
    beta= 2*(1-alpha**2)
    diagoprincipale=beta*np.ones(lvecteur+1)
    A=(np.diag(diagoprincipale)+np.diag(diagoinf,-1)+np.diag(diagoinf,1))
    A[0,1]=0
    A[lvecteur,lvecteur-1]=0
    A[0,0]=0
    A[lvecteur,lvecteur]=0
    
    U=np.zeros((lvecteur+1,cvecteur+1))
    #Resolution du prblm
    for i in range(2,cvecteur+1):
        U[:,i] = A @ U[:,i-1] - U[:,i-2] + (dt/2) * DiracF[:,i-1]
        U[0,i]=0
        U[lvecteur,i]=0
    return U
    
def Force(a,L,duree,f,amplitude,dt,dx,t):
    lvecteur=int(L/dx)
    cvecteur=int(duree/dt)
    periode=1/f
    diracIndice=int(a/dx)
    # etablissement de la matrice contenant le sinus porte
    dirac=np.array([0 for i in range(lvecteur+1)])
    dirac[diracIndice]=1
    dirac=dirac.reshape(lvecteur+1,1)
    
    F=amplitude*np.sin(2* np.pi * t / periode) * (t <= periode)
    F=F.reshape(1,cvecteur+1)
    DiracF=np.dot(dirac,F)
    
    return DiracF

t = np.linspace(0, duree, cvecteur+1)
F=Force(x,Longueur,duree,f,amplitude,pastemps,pasespace,t)
U=Matrice_deplacement(c,duree,Longueur,F,pastemps,pasespace)


# # # Observation en un point quelconque 

# plt.figure()

plt.plot(t,U[int(exemplepoint//pasespace),:],"b",label=f'amplitude en x={exemplepoint}')
plt.xlabel("temps")
plt.ylabel("déplacement")
plt.title(f'Solution FDM en x={exemplepoint} ')
plt.legend()
plt.grid(True)
plt.show()
# # Visualisation de la solution exacte
# plt.figure(figsize=(10, 6))
# plt.title("Solution différence finies de l'équation de d'Alembert")
# plt.pcolormesh(t, X, U, shading='auto', cmap='jet')
# plt.colorbar(label="Amplitude")
# plt.xlabel("Temps (s)")
# plt.ylabel("Position (x)")
# plt.grid(True)
# plt.show()
# plt.figure()
# plt.plot(t,DiracF[diracIndice,:],"r",label="sinus porte")
# plt.xlabel("temps")
# plt.ylabel("dépmacement")
# plt.title("probleme 1D propagation onde")
# plt.legend()
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
# ani = animation.FuncAnimation(fig, animate, frames=range(0, cvecteur, 100), interval=10, blit=True, repeat=False)

# # Enregistrement de l'animation
# ani.save("propagation_ondes_FDM.mp4", writer="ffmpeg", fps=30)
# plt.show()

# fig, ax=plt.subplots(subplot_kw={"projection":"3d"})
# x,y=np.meshgrid(t,X)
# surf=ax.plot_surface(x,y,U)
# plt.show()
