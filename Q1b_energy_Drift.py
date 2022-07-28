#Part b 
#Energy drift of Verlet

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
xprev, yprev, zprev =  x_init(Npart, vx, vy ,vz, x, y, z, dt)
Es = []

for step in range(NSteps):
    sum2 = 0
    fx, fy, fz, en = forcener(Npart, x, y, z, L, rcut)
    for i in range(Npart):
        xx = 2*x[i] - xprev[i] + (dt**2)*fx[i]
        yy = 2*y[i] - yprev[i] + (dt**2)*fy[i]
        zz = 2*z[i] - zprev[i] + (dt**2)*fz[i]
        vx[i] = (xx - xprev[i])/(2*dt)
        vy[i] = (yy - yprev[i])/(2*dt)
        vz[i] = (zz - zprev[i])/(2*dt)
        sum2 = sum2 + vx[i]**2 + vy[i]**2 + vz[i]**2
        xprev[i] = x[i]
        yprev[i] = y[i]
        zprev[i] = z[i]
        x[i] = xx
        y[i] = yy
        z[i] = zz
    Temp = sum2/(3*Npart)
    Etot = (en+0.5*sum2)/Npart
    Es.append(Etot)
delE = np.mean(np.abs(Es[0]-Es))/Es[0]
delEs = [np.mean(np.abs(Es[0]-Es[:i+1]))/Es[0] for i in range(NSteps)]
print('Energy drift = '+str(delE))
plt.plot(delEs)
plt.xlabel('Number of Steps')
plt.ylabel('Energy Drift')
plt.title('Energy Drift with number of steps for Verlet Algorithm')
plt.show()