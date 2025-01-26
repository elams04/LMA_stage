import numpy as np
import matplotlib.pyplot as plt
import time
import math
def SimulationFDM (c,pastemps,pasespace,f) :
    alpha=c*pastemps/pasespace
    
    Longueur=30
    duree=50

    f=0.25
    periode=1/f
    x=Longueur/2
    amplitude=100
    diracIndice=int(x//pasespace)
    lvecteur=int(Longueur/pasespace)
    cvecteur=int(duree/pastemps)

    # conditions de Neuman u0 deplacement à t=0  u1 etant la vitesse à t=0
    def u0 (x):
        return 0
    def u1 (x):
        return 0
    # Conditions de bords
    a=0
    b=0

    #etablissement de la matrice A
    gamma=alpha**2
    diagoinf=gamma*np.ones(lvecteur-1)
    beta=beta = 2*(1-alpha**2)
    diagoprincipale=beta*np.ones(lvecteur)
    A=(np.diag(diagoprincipale)+np.diag(diagoinf,-1)+np.diag(diagoinf,1))
    A[0,1]=0
    A[lvecteur-1,lvecteur-2]=0
    A[0,0]=a
    A[lvecteur-1,lvecteur-1]=b





    # etablissement de la matrice contenant le sinus porte
    dirac=np.array([0 for i in range(lvecteur)])
    dirac[diracIndice]=1
    dirac=dirac.reshape(lvecteur,1)
    t = np.linspace(0, duree, cvecteur)
    F=amplitude*np.sin(2* np.pi * t / periode) * (t <= periode)
    F=F.reshape(1,cvecteur)
    DiracF=np.dot(dirac,F)

    # à t=0 Vecteurs deplacement  U0 et vecteur vitesse V1
    U0=np.array([u0(i*pasespace) for i in range (lvecteur)])
    V1=np.array([u1(i*pasespace) for i in range(lvecteur)])


    U=np.zeros((lvecteur,cvecteur))
    U[:,0]=U0
    U[:,1]=U0+pastemps*V1
    X=[i*pasespace for i in range(lvecteur)]

    #Resolution du prblm
    for i in range(2,cvecteur):
        U[:,i]=np.dot(A,U[:,i-1])-U[:,i-2]+(pastemps**2)*DiracF[:,i-1]
    norme_L2 = np.linalg.norm(U,ord=2)
    return math.log(norme_L2, 10)
    

def SimulationNewmark(c,dt,dx,f):
    
    Longueur=30
    duree=18

    f=0.5
    periode=1/f
    x=Longueur/2
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
    norme_L2 = np.linalg.norm(U,ord=2)
    return math.log(norme_L2, 10)




t=[i*0.001 for i in range(1,120,5)]
alpha=[i/0.1 for i in t]
MAXFDM=[]
MAXNewmark=[]
TempsFDM=[]
TempsNewmark=[]
for dt in t:
    t1=time.perf_counter()
    MAXFDM.append(SimulationFDM(1,dt,0.1,0.5))
    t2=time.perf_counter()
    TempsFDM.append(t2-t1)
for dt in t:
    t1=time.perf_counter()
    MAXNewmark.append(SimulationNewmark(1,dt,0.1,0.5))
    t2=time.perf_counter()
    TempsNewmark.append(t2-t1)

plt.figure()
plt.plot(alpha,MAXFDM,'x',label='FDM')
plt.plot(alpha,MAXNewmark,'x',label='Newmark')
plt.xlabel("alpha")
plt.ylabel("log(L2)")

# plt.legend()
# plt.show()
# plt.figure()
# plt.plot(t,TempsFDM,label="FDM")
# plt.plot(t,TempsNewmark,label="Newmark")
# plt.xlabel("secondes")
# plt.xlabel("secondes")
# plt.legend()
# plt.figure()