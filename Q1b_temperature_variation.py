#Variation with Temp

from Q1a import *

dt = 0.00001
rhos = 1
Temps = [ 1000, 100, 10, 1, 0.1]
NSteps = 500
delEs = []
for Temp in Temps:
    delEs.append(sim_verlet(dt, rho, Temp, NSteps))
plt.plot(Temps, delEs)
plt.xlabel('Temperature')
plt.ylabel('Energy Drift')
plt.title('Energy Drift variation with Temperature')
plt.show()