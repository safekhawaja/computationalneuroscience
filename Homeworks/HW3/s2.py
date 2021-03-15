from mat4py import loadmat
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# a) Plot data from MATLAB file. We expect that the baseline voltages to be around -65 mV.

mats = loadmat('orientation_tuning_data.mat')

stimuli = mats['Stimuli']
info = mats['information']
vm = mats['Vm']

for key in info.keys():
    print(key, ":", info[key])

plt.plot(range(1477650), vm, color='#444444', linestyle='--', label='Vm')
plt.xlabel('Time (s)')
plt.ylabel('V (mV)')
plt.title('Regular Spiking Neurons')
plt.legend()
plt.show()

# b) Count the number of times that the voltage crosses the threshold to find firing rate.

counter = 0

t = 1
dt = 1

vthresh = -60

while t < 1477650:
    if vm[t] > vthresh:
        counter = counter + 1

print("In the data we have ", counter, " spikes.")

firing_rate = counter / 1477650

# c) c) Use # of spikes to plot firing rate against angle orientation (likely to find difference of 180 degrees).
# Plot a histogram for Spike triggered averages in response to different angles.

print("While I was unable to get this working, I would imagine that 0 and 180 degrees (given that they are horizontal, are the neurons' ideal angles given firing rates.")
# d) Use MATLAB function to fit Gaussian and find maximum (can also find mean and stdev to help fit).
# https://blog.dominodatalab.com/fitting-gaussian-process-models-python/

n = 1477650
mean = firing_rate / counter
sigma = sum(y * (angle - mean) ** 2) / n

x = np.random.randn(100000) * 50 + 75
x = np.round(x / 10) * 10
x = x[x >= 20]

yhist, xhist = np.histogram(x, bins=np.arange(4096))

xh = np.where(yhist > 0)[0]
yh = yhist[xh]


def gaussian(x, a, mean, sigma):
    return a * np.exp(-((x - mean) ** 2 / (2 * sigma ** 2)))


popt, pcov = curve_fit(gaussian, xh, yh, [10000, 100, 10])

plt.plot(yhist)
i = np.linspace(0, 300, 301)
plt.plot(i, gaussian(i, *popt))
plt.xlim(0, 300)
