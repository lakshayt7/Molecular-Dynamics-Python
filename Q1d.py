#Part d 

from Q1a import *

dt = 0.0001
Npart = 100
NSteps = 1000

sigma = 1
epsilon = 1
rcut = 2.5
rho = 1
Vol = Npart/rho
L = Vol**(1/3)
Temp = 10

x,y,z = Lattice(Npart, L)
vx, vy, vz = velocity_init(Npart, Temp)
vxhalf, vyhalf, vzhalf = velocity_init(Npart, Temp)
Es = []
fx, fy, fz, en = forcener(Npart, x, y, z, L, rcut)
for step in range(NSteps):
    sum2 = 0
    for i in range(Npart):
        x[i] = x[i] + vx[i]*dt + (dt**2)*(fx[i]/2)
        y[i] = y[i] + vy[i]*dt + (dt**2)*(fy[i]/2)
        z[i] = z[i] + vz[i]*dt + (dt**2)*(fz[i]/2)        
        vxhalf[i] = vx[i] + fx[i]*dt/2
        vyhalf[i] = vy[i] + fy[i]*dt/2
        vzhalf[i] = vz[i] + fz[i]*dt/2
    fx, fy, fz, en = forcener(Npart, x, y, z, L, rcut)
    for i in range(Npart):
        vx[i] = vx[i] + fx[i]*dt/2
        vy[i] = vy[i] + fy[i]*dt/2
        vz[i] = vz[i] + fz[i]*dt/2
        sum2 = sum2 + vx[i]**2 + vy[i]**2 + vz[i]**2
    Temp = sum2/(3*Npart)
    Etot = (en+0.5*sum2)/Npart
    Es.append(Etot)
delE = np.mean(np.abs((Es[0]-Es)/Es[0]))
delEs = [np.mean(np.abs((Es[0]-Es[:i+1])/Es[0])) for i in range(NSteps)]
print('Energy drift = '+str(delE))
plt.plot(delEs)
plt.xlabel('Number of Steps')
plt.ylabel('Energy Drift')
plt.title('Energy Drift with number of steps for Verlet Algorithm')
plt.show()