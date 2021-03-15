# -*- coding: utf-8 -*-
"""phys

# Homework 3

### Saif Khawaja 87474921

### Section 1

1a) u is the recovery variable that models inactivation of the sodium channels and opening of potassium channels
(in dv/dt, u provides the negative feedback and there is a sudden increase in u after a spike, making the same change in V require a greater injected current I).

b) c: the initial value that V returns to after the spike; if d increases, then u resets to greater value after a spike
and to increase V by the same amount we need to have a greater injection current I. This simulates the neuron’s refractory period.

c) du/dt is useful modelling the recovery period of the membrane, where a is the recovery timescale and b is the sensitivity of recovery to the membrane potential.
Increasing b

a/b change the dynamics of the membrane recovery variable.  What effect does parameter b have on the dynamics ofu? Think through the long term behavior of u if the membrane voltage is clamped to a single value. What effect does the parameter a have on the dynamics;  in particular, how does increasing or decreasingthe value change the dynamics?
"""
import matplotlib.pyplot as plt

# d) Regular Spiking (RS) neuron. Simulate the Izhikevich model with the parameters a = 0.02,b = 0.2,c = −65mV,d = 8.
# Feed the neuron a step current, I = 0 when t < 0 and I = 10 when t > 0.
# You should find a temporally regular spiking pattern not unlike the Integrate-and-Fire model from the last problem set.

a = 0.02
b = 0.2
c = -65
d = 8

I = {}
V = {}

t = -10
dt = 0.0001
tf = 250

V[t] = -65
V_t = -65
u = b * V_t

while t < tf:
    if t < 0:
        I[t] = 0
        du = a * (b * V_t - u) * dt
        dv = (0.04 * V_t * V_t + 5 * V_t + 140 - u) * dt

    else:
        I[t] = 10
        du = a * (b * V_t - u) * dt
        dv = (0.04 * V_t * V_t + 5 * V_t + 140 - u + 10) * dt

    V_t = V_t + dv
    u = u + du

    if V_t > 30:
        V_t = c
        u = u + d

    V[t] = V_t

    t = t + dt

print(V.keys, ":", V.values())

plt.plot(list(I.keys()), list(I.values()), color='#000000', linestyle='-', label='I(t)')
plt.plot(list(V.keys()), list(V.values()), color='#444444', linestyle='-', label='V(t)')
plt.xlabel('Time (s)')
plt.ylabel('I/V (A, mV)')
plt.title('Regular Spiking Neurons')
plt.legend()
plt.show()

# e) Fast Spiking (FS) neuron. Simulate the Izhikevich model with the parameters a = 0.1, b = 0.2, c =−65mV, d =  2.
# Feed  the  neuron  a  step  current, I =  0  when t < 0  and I = 10  when t > 0.
# You should find a temporally regular spiking pattern like in (a) but the firing rate should be higher for this fast-spiking neuron.
# Explain qualitatively how the larger value of a and the smaller value of d contribute to the fast spiking.
# Hint: Start with the RS model and change just the parametera, and then just d, and see what these single changes doto the firing pattern seen in part (a).
# Then see what happens if you change both parameters together.

a = 0.1
b = 0.2
c = -65
d = 2

I = {}
V = {}

t = -10
dt = 0.0001
tf = 250

V[t] = -65
V_t = -65
u = b * V_t

while t < tf:
    if t < 0:
        I[t] = 0
        du = a * (b * V_t - u) * dt
        dv = (0.04 * V_t * V_t + 5 * V_t + 140 - u) * dt

    else:
        I[t] = 10
        du = a * (b * V_t - u) * dt
        dv = (0.04 * V_t * V_t + 5 * V_t + 140 - u + 10) * dt

    V_t = V_t + dv
    u = u + du

    if V_t > 30:
        V_t = c
        u = u + d

    V[t] = V_t

    t = t + dt

print(V.keys, ":", V.values())

plt.plot(list(I.keys()), list(I.values()), color='#000000', linestyle='-', label='I(t)')
plt.plot(list(V.keys()), list(V.values()), color='#444444', linestyle='-', label='V(t)')
plt.xlabel('Time (s)')
plt.ylabel('I/V (A, mV)')
plt.title('Fast Spiking Neurons')
plt.legend()
plt.show()

# f) Chattering (CH) neuron. Simulate the Izhikevich model with the parameters a = 0.02, b = 0.2, c = −50mV,d=2.
# Feed the neuron a step current,I=0whent<0andI=10whent>0. You should see a response pattern that is very different from the RS and FS neurons simulated above.
# Explain qualitatively how the differences in c and d as compared to the RS neuron lead to the pattern that you see.

a = 0.02
b = 0.2
c = -50
d = 2

I = {}
V = {}

t = -10
dt = 0.0001
tf = 250

V[t] = -65
V_t = -65
u = b * V_t

while t < tf:
    if t < 0:
        I[t] = 0
        du = a * (b * V_t - u) * dt
        dv = (0.04 * V_t * V_t + 5 * V_t + 140 - u) * dt

    else:
        I[t] = 10
        du = a * (b * V_t - u) * dt
        dv = (0.04 * V_t * V_t + 5 * V_t + 140 - u + 10) * dt

    V_t = V_t + dv
    u = u + du

    if V_t > 30:
        V_t = c
        u = u + d

    V[t] = V_t

    t = t + dt

print(V.keys, ":", V.values())

plt.plot(list(I.keys()), list(I.values()), color='#000000', linestyle='-', label='I(t)')
plt.plot(list(V.keys()), list(V.values()), color='#444444', linestyle='-', label='V(t)')
plt.xlabel('Time (s)')
plt.ylabel('I/V (A, mV)')
plt.title('Chattering Spiking Neurons')
plt.legend()
plt.show()