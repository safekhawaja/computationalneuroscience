import math
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt

# g) Square and Sine Wave Currents

t = np.linspace(0, 1, 1000, endpoint=True)

a = 0.02
b = 0.2
c = -65
d = 8

Isq = {}
Isi = {}
V1 = {}
V2 = {}

t = -10
dt = 0.0001
tf = 250

V1[t] = -65
V2[t] = -65
V_t1 = -65
V_t2 = -65

u1 = b * V_t1
u2 = b * V_t2

while t < tf:
    if t < 0:
        Isquare = signal.square(2 * np.pi * 5 * t) * 2
        Isq[t] = Isquare
        Isine = math.sin(t / tf * 2 * math.pi) * 5
        Isi[t] = Isine
        du1 = a * (b * V_t1 - u1) * dt
        du2 = a * (b * V_t1 - u2) * dt
        dv1 = (0.04 * V_t1 * V_t1 + 5 * V_t1 + 140 - u1 + Isquare) * dt
        dv2 = (0.04 * V_t2 * V_t2 + 5 * V_t2 + 140 - u2 + Isine) * dt

    else:
        Isq[t] = signal.square(2 * np.pi * 5 * t) * 2
        Isi[t] = math.sin(t / tf * 2 * math.pi) * 5
        du1 = a * (b * V_t1 - u1) * dt
        du2 = a * (b * V_t1 - u2) * dt
        dv1 = (0.04 * V_t1 * V_t1 + 5 * V_t1 + 140 - u1 + Isquare) * dt
        dv2 = (0.04 * V_t2 * V_t2 + 5 * V_t2 + 140 - u2 + Isine) * dt

    V_t1 = V_t1 + dv1
    V_t2 = V_t2 + dv2

    u1 = u1 + du1
    u2 = u2 + du2

    if V_t1 > 30:
        V_t1 = c
        u1 = u1 + d

    if V_t2 > 30:
        V_t2 = c
        u2 = u2 + d

    V1[t] = V_t1
    V2[t] = V_t2

    t = t + dt

plt.plot(list(Isq.keys()), list(Isq.values()), color='#000000', linestyle='-', label='Square I(t)')
plt.plot(list(V1.keys()), list(V1.values()), color='#444444', linestyle='-', label='V(t)')
plt.xlabel('Time (s)')
plt.ylabel('I/V (A, mV)')
plt.title('Regular Spiking Neurons')
plt.legend()
plt.show()
plt.plot(list(Isi.keys()), list(Isi.values()), color='#000000', linestyle='-', label='Sine I(t)')
plt.plot(list(V2.keys()), list(V2.values()), color='#444444', linestyle='-', label='V(t)')
plt.xlabel('Time (s)')
plt.ylabel('I/V (A, mV)')
plt.title('Regular Spiking Neurons')
plt.legend()
plt.show()
