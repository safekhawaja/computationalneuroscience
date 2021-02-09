from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

'''
Try two other choices for Ie(t) (A periodic function for example).  
Add the graphs showing the current and voltage traces to your notes.
'''

# Collaborated with Cat G and Raveen K

# Sampling rate 1000 hz / second

t = np.linspace(0, 1, 1000, endpoint=True)

# Set constants

R_membrane = 1 * 10 ** 7-
V_threshold = -50 * (10 ** -3)
tau = 0.01
dt = 1 * (10**-5)
EL = -0.065
dv = 0

# Initialize dictionaries to place calculated values

V = {}
I_e = {}

t = 0
V[t] = EL
t_f = 1

# Obtain potentials by scanning across time range

while t < t_f:
    I_e_square = signal.square(2 * np.pi * 5 * t) * 2 * (10 ** -9)
    I_e[t] = I_e_square
    dv = (dt / tau) * ((EL - V[t]) + (R_membrane * I_e_square))
    V_temp = V[t] + dv
    t = t + dt

    V[t] = V_temp

    if V[t] >= V_threshold:
        V[t] = EL

plt.plot(list(V.keys()), list(V.values()), color='#444444', linestyle='--', label='')
plt.xlabel('Time (s)')
plt.ylabel('Potential (mV)')
plt.title('Square Regression')
plt.show()

plt.plot(list(I_e.keys()), list(I_e.values()), color='#444444', linestyle='--', label='')
plt.xlabel('Time (s)')
plt.ylabel('Current (nA)')
plt.title('Square Currents')
plt.show()