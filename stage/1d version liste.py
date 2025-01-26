# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 09:33:02 2024

@author: Utilisateur
"""
import numpy as np
import matplotlib.pyplot as plt
f=0.5
c=1
periode=1/f
dt=0.001
dx=0.1
Longueur=10
duree=30
x=1
xdirac=5
amplitude=100
N=int(duree//dt)
L=int(Longueur//dx)
alpha=c*dt/dx
print("alpha=",alpha)
beta=2*(1-alpha**2)
gamma=alpha**2
exemple=int(x//dx)
diracIndice=int(xdirac//dx)
point=[0,0]
t = np.linspace(0, duree, N)
F=amplitude*np.sin(2* np.pi * t / periode) * (t <= periode)
Un=[0 for i in range(L)]
Un_=[0 for i in range(L)]
U=[0 for i in range(L)]
for n in range(2,N):
    for l in range (1,L-1):
        if l==diracIndice:
            U[l]=Un[l]*beta-Un_[l]+gamma*(Un[l+1]+Un[l-1])+(dt*dt*F[n])
        else:
            U[l]=Un[l]*beta-Un_[l]+gamma*(Un[l+1]+Un[l-1])
        
    U[0]=0
    U[L-1]=0
    Un_=Un[:]
    Un=U[:]
    point.append(U[exemple])
# plt.figure()
# plt.plot(t,F)
# plt.show()
plt.figure()
plt.plot(t,point)
plt.grid(True)
plt.show()
    
