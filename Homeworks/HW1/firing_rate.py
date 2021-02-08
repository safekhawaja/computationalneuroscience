import matplotlib.pyplot as plt
import numpy as np

'''
Try different constant values for Ie and produce a graph showing how the firing rate changes withIe.
(The firing rate is defined as the number of threshold crossings, i.e., action potentials, per second.)  
You can either do this by hand or write a MATLAB program doing it for you.  
What is approximately the minimum value for Ie after which the cell starts to produce action potentials at all?
'''

# Collaborated with Cat G and Raveen K

# Set constants

I_e = 2 * (10 ** -9)
R_membrane = 1 * 10 ** 7
V_threshold = -50 * (10 ** -3)
tau = 0.01
dt = 1 * (10 ** -5)
EL = -0.065
dv = 0

currents = np.arange(1, 6, 0.5)
currents = currents * 10 ** (-9)

fire_rate = list()

dt = 1 * 10 ** (-6)
t_final = 1

V = {}

for i in currents:
    counter = 0

    t = 0
    V[t] = EL

    while t < t_final:
        dV = (1 / tau) * (EL - V[t] + R_membrane * i) * dt
        V_new = V[t] + dV

        t = t + dt

        if V_new > V_threshold:
            V[t] = EL
            counter += 1
        else:
            V[t] = V_new

    fire_rate.append(counter)

plt.plot(currents, fire_rate, color='#444444', linestyle='--', label='')
plt.xlabel('Current)')
plt.ylabel('Fire Rate')
plt.show()

print("Approximate value for action potentials is is 1.5 nA")
