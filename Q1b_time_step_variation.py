#Variation with time step
from Q1a import *

dts = [10**(-i-4) for i in range(6)]
Npart = 100
rho = 1
Temp = 1
NSteps = 100
delEs = []
for dt in dts:
    delEs.append(sim_verlet(dt, rho, Temp, NSteps))
plt.plot(np.log10(dts), np.log10(delEs))
plt.xlabel('log Time Step')
plt.ylabel('log Energy Drift')
plt.title('Energy Drift variation with Time Step on log log scale')
plt.show()