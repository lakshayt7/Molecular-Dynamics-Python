import math
import numpy as np
import matplotlib.pyplot as plt

sigma = 1
epsilon = 1
rcut = 2.5
rho = 1
Npart = 100
Vol = Npart/rho
L = Vol**(1/3)
Temp = 1

def Lattice(Npart, L):
    X = []
    Y = []
    Z = []
    
    N = int((Npart)**(1/3))+1
    if N == 0:
        N = 1     
    Itel = 0
    Del = L/float(N)
    Dx = -Del
    for I in range(0, N, 1):
        Dx = Dx + Del
        Dy = -Del
        for J in range(0, N, 1):
            Dy = Dy + Del
            Dz = -Del
            for K in range(0, N, 1):
                Dz = Dz + Del
                if Itel < Npart:
                    Itel = Itel + 1
                    X.append(Dx)
                    Y.append(Dy)
                    Z.append(Dz)              
    return X, Y, Z 

def velocity_init(Npart, Temp):
    vx = np.random.uniform(size=Npart) - 0.5
    vy = np.random.uniform(size=Npart) - 0.5
    vz = np.random.uniform(size=Npart) - 0.5
    sumx = np.mean(vx)
    sumy = np.mean(vy)
    sumz = np.mean(vz)
    sum2 = np.mean(np.multiply(vx, vx)+np.multiply(vy, vy)+np.multiply(vz, vz))
    fs = math.sqrt(3*Temp/sum2)
    vx = (vx-sumx)*fs
    vy = (vy-sumy)*fs
    vz = (vz-sumz)*fs
    return vx,vy,vz

def x_init(Npart, vx, vy ,vz, x, y, z, dt):
    return x - vx*dt, y - vy*dt, z - vz*dt

def forcener(Npart, x, y, z, L, rcut):
    en = 0
    fx = np.zeros((Npart))
    fy = np.zeros((Npart))
    fz = np.zeros((Npart))
    rci = 1/rcut
    rci6 = rci**6
    ecut = 4*(rci6**2 - rci6)
    for i in range(Npart-1):
        for j in range(i+1, Npart):
            xr = x[i] - x[j]
            yr = y[i] - y[j]
            zr = z[i] - z[j]
            xr = xr - L*(round(xr/L))
            yr = yr - L*(round(yr/L))
            zr = zr - L*(round(zr/L))
            r2 = (xr**2)+(yr**2)+(zr**2)
            if r2<(rcut**2):
                r2i = 1/r2
                r6i = r2i**3
                ff = 48*r2i*r6i*(r6i - 0.5) 
                fx[i] = fx[i] + ff*xr
                fy[i] = fy[i] + ff*yr
                fz[i] = fz[i] + ff*zr
                en = en + 4*r6i*(r6i - 1) - ecut
    return fx, fy, fz, en
#part b variation with dt 
def sim_verlet(dt, rho, Temp, NSteps = 1000):
    sigma = 1
    epsilon = 1
    rcut = 2.5
    rho = 1
    Vol = Npart/rho
    L = Vol**(1/3)
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
    return delE

if __name__ == "__main__":
    dt = float(input("Enter time step dt : "))
    Npart = int(input("Enter number of particles : "))
    NSteps = int(input("Enter number of steps : "))
    x,y,z = Lattice(Npart, L)
    vx, vy, vz = velocity_init(Npart, Temp)
    xprev, yprev, zprev =  x_init(Npart, vx, vy ,vz, x, y, z, dt)
    Es = []

    for step in range(NSteps):
        sum2 = 0
        try: 
            fx, fy, fz, en = forcener(Npart, x, y, z, L, rcut)
        except:
            print("Simulation diverged to nan (infinity) values because time step is too high! Reduce time step.")
            raise
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
    print('Energy = ' + str(np.mean(Es)))
    print('Temperature = ' + str(Temp)) 