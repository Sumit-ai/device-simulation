##############################################################################################################
'''                       PART 1 Constents initialization             '''
#############################################################################################################

from math import pi
import math 
## Constants for intrinsic carrier conc.
#pi=3.14
Kb=8.617*10**(-5) #Boltmann constant in eV/K
h=4.13*10**(-15)  #Planck constant in eV/s
mn=1.1*9.1*10**(-31) #Effective mass of electron in Kg
mp=0.59*9.1**(-31) #Effective mass of holes in Kg
Eg=1.1 #Bandgap of intrinsic Silicon
######
## Constants for Boron Diffusion ### 
Dixo=0.037 #in cm2/sec  
Dipo=0.72 #in cm2/sec
Dixe=3.46 #in eV
Dipe=3.46 #in eV
######
T = 1323 # 1050 C in Kelvin 
ni=((2*pi*Kb*T/(h**2))**1.5)*((mn*mp)**(3/4))*math.exp(-Eg/(2*Kb*T))
Dix=Dixo*math.exp(-Dixe/(Kb*T))
Dip=Dipo*math.exp(-Dipe/(Kb*T))
D=Dix+Dip*(ni/(10**18))

print('value of D at temperature: ', T, ' is ', D)

##############################################################################################################
'''                       PART 2 backward euler              '''
#############################################################################################################


import numpy as np
import matplotlib.pyplot as plt 
plt.grid()
# intialize all the constents here
#D = 7.36754309270159e-15
dx = 1.35*10**(-7) # make sure agjust this value to get the value of alpha as less or equal to 0.45

L = 24;S= 100;alpha=D/dx**2;T=60;

x = np.linspace(0, L, S+1)   # mesh points in space
T = 60
t = np.linspace(0, T, T+1)   # mesh points in time
Sp   = np.zeros(S+1)

a = 10**18*np.ones(40)
b = 10**15*np.ones(61)
space = np.concatenate((a,b))

for n in range(1, T):
    '''Compute Sp at inner mesh points'''
    for i in range(1, S):
        Sp[i] = space[i] + alpha*(space[i-1] - 2*space[i] + space[i+1])
    # Insert boundary conditions
    Sp[0] = 10**18;  Sp[S] = 10**15
    # Switch variables before next step
    space, Sp = Sp, space
    #plt.clf()
    plt.subplot(211)
    
    plt.plot(x, Sp)
    
    SRF = "Diffusion of Boron in silicon after " + str(n+1) + ' minutes'
    plt.title(SRF)
    plt.ylabel("Concetration")
    plt.xlabel("length (in micrometer)")
    plt.subplot(212)
    plt.plot(x, np.log(Sp))
    
    PRF = "Diffusion of Boron in silicon after " + str(n+1) + " minutes (log plot)"  
    plt.title(PRF)
    plt.ylabel("Concetration")
    plt.xlabel("length (in micrometer)")
    plt.tight_layout()
    plt.pause(0.2)
    #plt.plot(u_1)

