%% 1a
load('contrast_response.mat')
sampleRate = 2000;
totalInterval = 20;
numSamples = sampleRate*totalInterval;

spikeTimesSelected = spikeTimes{9};
mask = spikeTimesSelected<=numSamples;
actionPotentials = mask(mask~=0);
times = spikeTimesSelected(mask)/2000;

figure
stem(times, actionPotentials);
xlabel('Time (s)');
ylabel('Action Potential (1 for spike, 0 for no spike)');
title('Spikes for first 20 seconds of trial with highest contrast')

%% 1b
[numSpikes, sampleNumbers] = hist(spikeTimesSelected(mask), 200);
rate = numSpikes/0.1;
times = sampleNumbers/2000;

figure
plot(times,rate)
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Time-dependent firing rate for first 20 seconds of trial with highest contrast')
grid on

%% 1c
x = -5000:1:5000;
sigma = 70;
gauss = gaussmf(x, [sigma 0]);
gauss = 2000*gauss/(sum(gauss));

for i = 1:9
    contrast = contrasts(i);
    spikeTimesSelected = spikeTimes{i};
    spikeArray = zeros(1, numel(stimulus));
    spikeArray(spikeTimesSelected) = 1;
    convnFiring{i} = conv(spikeArray(1:numel(stimulus)), gauss, 'same');
    time = 1:numel(stimulus);
  
    subplot(3, 3, i);
    plot(time/2000,convnFiring{i})
    xlabel('Time (s)');
    ylabel('Firing Rate (spikes/s)');
    title(['Contrast = ', num2str(contrast)]);
    grid on
end

%% 1d
firingRate = zeros(1, 9);

for i = 1:9
    numSpikes = numel(spikeTimes{1,i});
    firingRate(i) = numSpikes/100;
end

plot(contrasts, firingRate)
xlabel('Contrast');
ylabel('Firing Rate (spikes/s)');
title('Average Firing Rate as a function of Contrast')
grid on

%% 1e
tau = 0:1:500;
y = zeros(numel(tau),9);

for i = 1:9
    contrast = contrasts(i);
    stimVector = stimulus * contrast;
    spikeTimesSelected = spikeTimes{1,i};
    for t=1:numel(tau)
        spikeTimesSelectedNew = spikeTimesSelected - tau(t);
        spikeTimesSelectedNew = spikeTimesSelectedNew(spikeTimesSelectedNew>=1);
        stimValues = stimVector(spikeTimesSelectedNew); 
        stimValuesSummed = sum(stimValues);
        y(t,i) = stimValuesSummed/numel(spikeTimesSelected);
    end
end

for i = 1:9
    y(:,i) = y(:,i)./(norm(y(:,i)));  
end

figure
for i = 1:9
    contrast = contrasts(i);
    subplot(3,3,i);
    plot(tau/2, y(:,i));
    xlabel('\tau (ms)');
    ylabel('Stimulus Value');
    title(['Contrast = ', num2str(contrast)]);
    xlim([0,250])
    set(gca, 'XDir','reverse')
    grid on
end

%% 2a
for i = 1:9
    contrast = contrasts(i);
    stimVector = stimulus * contrast;
    STA = y(:,i);
    triggerIntensity{i} = conv(stimVector, STA, 'full');
end

for i = 1:9
    contrast = contrasts(i);
    subplot(3,3,i);
    times = (1:numel(triggerIntensity{i}))/numel(triggerIntensity{i})*20;
    plot(times, triggerIntensity{i});
    xlabel('Time (s)');
    ylabel('Trigger Feature Intensity');
    title(['Contrast = ', num2str(contrast)]);
    grid on
end

%% 2b
contrastTrialNum = 9; % change this value to affect which contrast trial will be plotted, need to run this once for every trial
sigma = 20;
time = 1:numel(stimulus); 

x = -5000:.5:5000;
gauss = gaussmf(x, [sigma 0]);
gauss = 2000*gauss/(sum(gauss));
spikeTimesSelected = spikeTimes{contrastTrialNum};
spikeArray = zeros(1, numel(stimulus));
spikeArray(spikeTimesSelected) = 1;

fireFrom1c = conv(spikeArray, gauss, 'same');
fireFrom1c = fireFrom1c(1:numel(stimulus));           

triggerFrom2a = triggerIntensity{contrastTrialNum};
triggerFrom2a = triggerFrom2a(1:numel(stimulus));
triggerFrom2a = triggerFrom2a*(max(fireFrom1c)/max(triggerFrom2a));

triggerStore{contrastTrialNum} = triggerFrom2a;
fireStore{contrastTrialNum} = fireFrom1c;

triggerFrom2a(triggerFrom2a<0)=0;

figure
plot(time/2000,fireFrom1c)
hold on
plot(time/2000,triggerFrom2a)
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title(['Trial #',num2str(contrastTrialNum)]);
xlim([10,15])
legend('True Firing Rate','Trigger Feature Intensity')
grid on

%% 2c
for i = 1:9
    trigIntensity = triggerStore{i};
    firingRate = fireStore{i};
    
    subplot(3,3,i)
    scatter(trigIntensity(1:10:end), firingRate(1:10:end))
    xlabel('Trigger Feature Intensity');
    ylabel('Firing Rate (spikes/s)');
    title(['Trial #',num2str(i)]);
    grid on
end

%% 2d
for i = 1:9
    trigIntensity = transpose(triggerStore{i});
    firingRate = (fireStore{i});
    [N, edges] = histcounts(trigIntensity);

    for j = 1:numel(edges)-1
        edgeLower = edges(j);
        edgeUpper = edges(j+1);
        binnedIndices = find((trigIntensity>=edgeLower)&(trigIntensity<=edgeUpper));
        correspondingFiringRates = firingRate(binnedIndices);
        avgFiringRate(j) = mean(correspondingFiringRates);
    end
    
    midpt = (edges(2)-edges(1))/2;
    filterVals = edges-midpt;
    filterVals = filterVals(1:(numel(edges)-1));
    filterVals = filterVals(1:3:end);
    avgFiringRate = avgFiringRate(1:3:end);
    filterValsTrial{i} = filterVals;
    avgFiringRateTrial{i} = avgFiringRate;

    subplot(3,3,i)
    scatter(filterVals, avgFiringRate)
    xlabel('Trigger Feature Intensity');
    ylabel('Firing Rate (spikes/s)');
    title(['Trial #', num2str(i)]);
    grid on
end

%% 2e
for i = 1:9
    predictFire = interp1(filterValsTrial{i}, avgFiringRateTrial{i}, triggerStore{i});
    times = (1:numel(predictFire))/2000;
    
    subplot(3,3,i)
    plot(times, fireStore{i})
    hold on
    plot(times,predictFire)
    xlabel('Time (s)');
    ylabel('Firing Rate (spikes/s)');
    title(['Trial #',num2str(i)]);
    xlim([10,13])
    legend('Actual','Prediction')
    set(gcf,'position',[0,0,900,600])
    grid on
end

%% 3a
phi = 0;
k = 1;
sigma = 2/k;
sigmax = sigma;
sigmay = sigma;

figure
[x, y] = meshgrid(-5:.1:5);
z = (1/(2*pi*sigmax*sigmay))*exp(-(x.^2)/(2*(sigmax.^2))-(y.^2)/(2*(sigmay.^2))).*cos(k*x - phi);
surf(x,y,z)
xlabel('x (degrees)')
ylabel('y (degrees)')
title('Gabor Function as 2D Surface');

figure
image(z,'CDataMapping','scaled');
set(gca,'XTickLabel',-4:5)
set(gca,'yTickLabel',-4:5)
xlabel('x (degrees)')
ylabel('y (degrees)')
title('Gabor Function as Colorplot');

%% 3b
tau = 0:0.5:300;
alpha = 1/15;
vals = alpha*exp(-alpha*tau).*(((alpha*tau).^5)/(125)-((alpha*tau).^7)/5040);

plot(tau,vals)
set(gca,'XDir','reverse')
xlabel('\tau (ms)')
ylabel('Dt (Hz)')
title('Temporal Evolution of Spatial Receptive Field');

%% 3c
phi = 0;
k = 1;
sigma = 2/k;
sigmax = sigma;
sigmay = sigma;

tau = 50:50:300;
alpha = 1/15;
[x, y] = meshgrid(-5:.1:5);

for i=1:numel(tau)
     z = ((1/(2*pi*sigmax*sigmay))*exp(-(x.^2)/(2*(sigmax.^2))-(y.^2)/(2*(sigmay.^2))).*cos(k*x-phi)).*(alpha*exp(-alpha*tau(i)).*(((alpha*tau(i)).^5)/(125)-((alpha*tau(i)).^7)/5040));
     subplot(2,3,i);  
     contour(x,y,z)
     xlabel('x (degrees)')
     ylabel('y (degrees)')
     title(['\tau = ', num2str(tau(i)), 'ms']);
     colorbar('eastoutside')
end

%% 3d
t = 0;
phi = 0;
A = 1;
w=800;

x = 0:0.1:200;
y = 0:0.1:200;
K = 0.2;
theta = 0.5*pi;

z = A*cos((K*x*cos(theta))+(K*y*sin(theta))-phi).*cos(w*t);

figure
image(z,'CDataMapping','scaled');
title('Counterphase Grating for t = 0');
set(gca,  'xtick', [])
set(gca, 'ytick', [])

t = linspace(0, 2, 20);
theta = linspace(0, 2*pi, 20);

figure
for i=1:numel(t)
    z = A*cos((K*x*cos(theta(i)))+(K*y*sin(theta(i)))-phi).*cos(w*t(i));
    subplot(4,5,i);
    image(z,'CDataMapping','scaled');
    title(['t = ', num2str(t(i))]);
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
end

figure
for i=1:numel(t)
    z = A*cos((K*x*cos(theta(i)))+(K*y*sin(theta(i)))-phi).*cos(w*t(i));
    image(z,'CDataMapping','scaled');
    title(['Counterphase Grating for t = ', num2str(t(i))]);
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
    pause(0.3);
end

%% 3e
clear variables
k = 1;
sigma = 2/k;
sigmax = sigma;
sigmay = sigma;
A = 1;
phi = 0;
K =k;
theta = -1.5:.02:1.5;
for i=1:numel(theta)
    fun = @(x,y) ((1/(2*pi*sigmax*sigmay))*exp(-(x.^2)/(2*(sigmax.^2))-(y.^2)/(2*(sigmay.^2))).*cos(k*x-phi)).*(A*cos(K.*x.*cos(theta(i))+K.*y.*sin(theta(i))-phi));
    Ls_A(i) = quad2d(fun,-100,100,-100,100);
end
figure
plot(theta, Ls_A)
xlabel('\theta (radians)')
ylabel('L_{s}')
title('Orientation Response (Figure A)');

clear variables
theta = 0;
k = 1;
sigma = 2/k;
sigmax = sigma;
sigmay = sigma;
A = 1;
phi = 0;
K = 0:0.02:3;
for i=1:numel(K)
    fun = @(x,y) ((1/(2*pi*sigmax*sigmay))*exp(-(x.^2)/(2*(sigmax.^2))-(y.^2)/(2*(sigmay.^2))).*cos(k*x-phi)).*(A*cos(K(i).*x.*cos(theta)+K(i).*y.*sin(theta)-phi));
    Ls_B(i) = quad2d(fun,-100,100,-100,100);
end
figure
plot(K, Ls_B)
xlabel('K/k')
ylabel('L_{s}')
title('Spatial Frequency Response (Figure B)');

clear variables
theta = 0;
k = 1;
sigma = 2/k;
sigmax = sigma;
sigmay = sigma;
A = 1;
K = k;
phi = -2.5:0.02:2.5;
for i=1:numel(phi)
    fun = @(x,y) ((1/(2*pi*sigmax*sigmay))*exp(-(x.^2)/(2*(sigmax.^2))-(y.^2)/(2*(sigmay.^2))).*cos(k*x-phi(i))).*(A*cos(K.*x.*cos(theta)+K.*y.*sin(theta)-phi(i)));
    Ls_C(i) = quad2d(fun,-50,50,-50,50);
end
figure
plot(phi, Ls_C)
xlabel('\phi')
ylabel('L_{s}')
title('Spatial Phase Response (Figure C)');

%% 3f
clear variables
alpha = 1/15;
dt = 0.2;
tau = 0:dt:1000;
t = 0:dt:1000;
w = 0:.01*.02:2*pi * .02;

for i=1:numel(w)
    func1 = alpha*exp(-alpha*tau).*(((alpha*tau).^5)/(125)-((alpha*tau).^7)/5040);
    func2 = cos(w(i)*(t));
    convolution = dt*conv(func2,func1,'same');
    pks=max(convolution);
    amplitude(i) = (pks);
end

plot(w/(2*pi)*1000, amplitude)  
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Model Frequency Response');