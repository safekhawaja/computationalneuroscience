import matplotlib.pyplot as plt

'''
Now implement the Integrate-and-Fire model for constant input currents.  
Let Ie(t) be a positive constantin a reasonable range and set the parameters as described above.  
Choose dt to be small enough to get smooth results (this will depend on your choice of Ie(t)).  
Plot the voltage as a function of time including a number of voltage resets as the threshold potential was reached.  
(If you cannot get resets, your Ie is probably too small.)  
Comment your code properly and add it, along with all later code, to the printout of your homework!
'''

# Collaborated with Cat G and Raveen K

# Set constants

I_e = 2 * (10 ** -9)
R_membrane = 1 * 10 ** 7
V_threshold = -50 * (10 ** -3)
tau = 0.01
dt = 1 * (10**-5)
EL = -0.065
dv = 0

# Initialize dictionaries to place calculated values

V = {}

t = 0
V[t] = EL
t_f = 0.1

# Obtain potentials by scanning across time range

while t < t_f:
    dv = (dt / tau) * ((EL - V[t]) + (R_membrane * I_e))
    V_temp = V[t] + dv
    t = t + dt

    V[t] = V_temp

    if V[t] >= V_threshold:
        V[t] = EL

# Plot axes

plt.plot(list(V.keys()), list(V.values()), color='#444444', linestyle='--', label='')
plt.xlabel('Time (s)')
plt.ylabel('Potential (mV)')
plt.title('Test Regression')
plt.show()
