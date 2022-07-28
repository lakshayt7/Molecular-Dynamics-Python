#Variation with density

from Q1a import *

dt = 0.0001
rhos = [0.01, 1, 100]
Temp = 1
NSteps = 2000
delEs = []
for rho in rhos:
    delEs.append(sim_verlet(dt, rho, Temp, NSteps))
plt.plot(np.log10(rhos), delEs)
plt.xlabel('log Density')
plt.ylabel('Energy Drift')
plt.title('Energy Drift variation with log density')
plt.show()