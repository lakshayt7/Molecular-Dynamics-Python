#part c

from Q1a import *

def Pressure(Npart, x, y, z, L, rho, Temp, rcut):
    Pressure = 0
    Vol = L**3
    rci = 1/rcut
    rci3 = rci**3
    Beta = 1/Temp
    Ptail = (16*3.14/3)*rho*rho*( (2*rci3*rci3*rci/3) - rci3)
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
                Pressure = Pressure + 48*(r6i*r6i - r6i/2)/(3*Vol)
    Pressure = Pressure + Ptail +  rho/Beta
    return Pressure

def sim_pressure(rho, Temp, NSteps = 1000):
    dt = 0.00000001
    Vol = Npart/rho
    L = Vol**(1/3)
    x,y,z = Lattice(Npart, L)
    vx, vy, vz = velocity_init(Npart, Temp)
    xprev, yprev, zprev =  x_init(Npart, vx, vy ,vz, x, y, z, dt)
    Es = []
    Pressures = []
    for step in range(NSteps):
        sum2 = 0
        fx, fy, fz, en = forcener(Npart, x, y, z, L, rcut)
        Pres = Pressure(Npart, x, y, z, L, rho, Temp, rcut)
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
        Pressures.append(Pres)
    return np.mean(Pressures[:-int(NSteps/2)])
rhos = [0.1*(i+1) for i in range(10)]
Temps = [0.01, 1, 100]
for Temp in Temps:
    Pressures = []
    for rho in rhos:
        pres = sim_pressure(rho, Temp, 100)
        Pressures.append(pres)
    plt.plot(rhos, Pressures)
    plt.xlabel('Density')
    plt.ylabel('Pressure')   
    plt.title('Pressure versus Density plot for Temperature = ' + str(Temp))
    plt.show()