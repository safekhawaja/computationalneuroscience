import matplotlib.pyplot as plt

'''
Try different constant values forIeand produce a graph showing how the firing rate changes withIe.
(The firing rate is defined as the number of threshold crossings, i.e., action potentials, per second.)  
You can either do this by hand or write a MATLAB program doing it for you.  
What is approximately the minimum value for Ie after which the cell starts to produce action potentials at all?
'''

# Collaborated with Cat G and Raveen K

# Set constants

I_e = 2e-9;
V_threshold=-50e-3;
R_membrane=1*10^7;
tau=0.01;
dt=1e-5;

V=[];
p_vector=[];
V(1)=(membrane_potential / 1000);
currents = 1e-9:1e-10:4e-9;
for j = 1:length(currents) peaks = 0; for i = 2:length(range) dv = (dt/tau)((membrane_potential / 1000) - V(i-1) + (R_membranecurrents(j))); dv; V(i) = dv+V(i-1); if V(i) >= V_threshold peaks = peaks + 1; V(i)= membrane_potential / 1000; end
end
p_vector = [p_vector peaks];
end
figure

plt.plot(list(V.keys()), list(V.values()), color='#444444', linestyle='--', label='')
plt.xlabel('Time (s)')
plt.ylabel('Potential (mV)')
plt.title('Test Regression')
plt.show()