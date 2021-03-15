%% 1(a)
sampleRate = 2000;
totalInterval = 20;
numSamples = sampleRate*totalInterval;

spikeTimesSelected = spikeTimes{9};
mask = spikeTimesSelected<=numSamples;
actionPotentials = mask(mask~=0);
times = spikeTimesSelected(mask)/2000;

figure
plt = stem(times, actionPotentials);
set(plt, 'Marker', 'none')
ylim([0,2])
xlabel('Time (s)','FontSize',16);
ylabel('Action Potential (0 = No Spike,1 = Spike)','FontSize',16);
title('Action Potentials for First 20 Seconds','FontSize',16)
width=900;
height=600;
set(gcf,'position',[0,0,width,height])


%% Part 1(b)
[nSpikes_inBin,time] = hist(spikeTimesSelected(mask),200); %200 100ms time bins

rate = nSpikes_inBin/.1;

plot(time/2000,rate,'LineWidth',2)
ylim([-1,91])
xlabel('Time (s)','FontSize',16);
ylabel('Firing Rate (Counts/Sec)','FontSize',16);
title('Neuron Firing Rate for First 20 Seconds','FontSize',16)
width=900;
height=600;
set(gcf,'position',[0,0,width,height])

%% Part 1(c) 

x = 0:200000;
for t=1:9
    x = 0:200000;
    contrast = contrasts(t);
    spikeTimesSelected = spikeTimes{t}; %center gaussians on spikeTimesSelected
    
    for h=1:numel(spikeTimesSelected)
        center = spikeTimesSelected(h);
        sigma = 70;
     %   sigma = 20;
        gaussian = gaussmf(x,[sigma center]);
        gaussian = 2000*gaussian/(sum(gaussian));
        vector{h} = gaussian;
              
    end
    
    mysum = vector{1};
    for i = 2:numel(vector)
        mysum = mysum + vector{i};
        
    end
    subplot(3,3,t);  
    plot(x/2000,mysum)
    xlim([0,100]);
    xlabel('time (s)','FontSize',16);
    ylabel('Firing Rate (Hz)','FontSize',16);
    title(['contrast =',' ',num2str(contrast)]);
    width=900;
    height=600;
    set(gcf,'position',[0,0,width,height])
    
    clearvars -except t
    load('contrast_response.mat')
    
end
    
    
%% Part 1(c) attempted again but with a convolution USE THIS.... PLOTS ARE SAME CODE IS BETTER
x = -5000:1:5000;
sigma = 70;
gaussian = gaussmf(x,[sigma 0]);
gaussian = 2000*gaussian/(sum(gaussian));
for trial = 1:9
    contrast = contrasts(trial);
    spikeTimesSelected = spikeTimes{trial};
    spike_array = zeros(1,numel(stimulus));
    spike_array(spikeTimesSelected)=1;
    convolution_firing{trial} = conv(spike_array(1:numel(stimulus)),gaussian,'same'); %trial elemnt of array is convolution values. We will use this later!
    time = 1:numel(stimulus);
    
    
    subplot(3,3,trial);
    plot(time/2000,convolution_firing{trial})
    xlabel('time (s)','FontSize',16);
    ylabel('Firing Rate (Hz)','FontSize',16);
    title(['contrast =',' ',num2str(contrast)]);
    width=900;
    height=600;
    set(gcf,'position',[0,0,width,height])
    
end


   
%% Part 1(d)
firing_rate = zeros(1,9);
for i=1:9
    nSpikes = numel(spikeTimes{1,i});
    firing_rate(i) = nSpikes/(100); %avg rate in a trial = number spikes/time of trial
    
end
plot(contrasts,firing_rate,'LineWidth',2)
xlabel('Contrast','FontSize',16);
ylabel('Firing Rate (Counts/Sec)','FontSize',16);
title('Average Firing Rate for each Contrast Step','FontSize',16)
width=900;
height=600;
set(gcf,'position',[0,0,width,height])


%% Part 1(e)

tau = 0:1:500; %(tau values from 0 to 250ms  in 1ms  steps [using index values])
values = zeros(numel(tau),9);
for c=1:9          %c for each contrast (1-9)
    contrast = contrasts(c);
    stim = stimulus * contrast; %stimulus vector for ith contrast (i: 1->9)
    spikeTimesSelected = spikeTimes{1,c};
    
    for t=1:numel(tau)
        indices_spiketime_minus_tau = (spikeTimesSelected) - tau(t);
        %Next line removes all index values less than 1 (out of range)
        indices_spiketime_minus_tau = indices_spiketime_minus_tau(indices_spiketime_minus_tau>=1);
        
        stim_values = stim(indices_spiketime_minus_tau); %Values of stimulus tau seconds before spike
        summed = sum(stim_values);  %Sum all stimulus values
        %Now, store summed value in matrix with tau rows and 9 columns (tau
        %in rows, trials in columns)
        values(t,c) = summed/numel(spikeTimesSelected); %columns have trials, rows have tau values of STA. Each row is a different tau.
    end
end

%Now we normalize each column of the values matrix...
for d=1:9
    values(:,d)=values(:,d)./(norm(values(:,d)));  
end



figure
for p=1:9
    contrast = contrasts(p);
    subplot(3,3,p);
   
    plot(tau/2,values(:,p),'LineWidth',2);
    title(['contrast =',' ',num2str(contrast)]);
    xlabel('\tau (ms)','FontSize',16);
    ylabel('Stimulus Value','FontSize',16);

    width=900;
    height=600;
    set(gcf,'position',[0,0,width,height])
    set(gca, 'XDir','reverse')  % reverse x dir (tau will decrease from left to right)
    ylim([-0.15,0.18])
end    


%% Please see part 1(f) in PDF document


%% Part 2(a)
for c=1:9          %c for each contrast (1-9)
    contrast = contrasts(c);
    stim = stimulus * contrast; %stimulus vector for ith contrast (i: 1->9)
    STA = values(:,c);
    convolution = conv(STA,stim,'full');
    start_index = 1;
    final_index = 1000;
    time = (start_index:1:final_index)/2000;
    
    
    subplot(3,3,c);
   
    plot(time,convolution(start_index:1:final_index),'LineWidth',2);
    xt = get(gca, 'XTick');                               
    %set(gca, 'XTick', xt, 'XTickLabel', xt*1000)
    title(['contrast =',' ',num2str(contrast)]);
    xlabel('time (ms)','FontSize',16);
    ylabel('Est. Firing Rate (Hz)','FontSize',16);
    width=900;
    height=600;
    set(gcf,'position',[0,0,width,height])
end

subplot(3,3,1)
ylabel('Est. Firing Rate (Hz)','FontSize',16);
subplot(3,3,4)
ylabel('Est. Firing Rate (Hz)','FontSize',16);
subplot(3,3,7)
ylabel('Est. Firing Rate (Hz)','FontSize',16);

subplot(3,3,7)
xlabel('time (ms)','FontSize',16);
subplot(3,3,8)
xlabel('time (ms)','FontSize',16);
subplot(3,3,9)
xlabel('time (ms)','FontSize',16);


%% Part 2(a) another attempt

for trial = 1:9
    contrast = contrasts(trial);
    stim = stimulus * contrast;
    STA = ((values(:,trial)));
    trigger_intensity{trial} = conv(stim,(STA),'full');
    
end


%Now a quick setup for part 2(b)...
filter_output = {};
adjusted_smooth_fire = {};

L = trigger_intensity{9};
R = convolution_firing{9};
% 
scatter(L(1:10:numel(stimulus)),R(1:10:end))

%% Part 2(b)
%First, set all values < 0 from plots in 2(a) (trigger feature intensity)
%to 0. This will be accomplished by rerunning the above code with slight
%modifications.

for c=9        %c for each contrast (1-9)
    contrast = contrasts(c);
    stim = stimulus * contrast; %stimulus vector for ith contrast (i: 1->9)
    STA = values(:,c);
    numel(STA)
    convolution = conv(STA,stim,'full');
    convolution = convolution(convolution>0); %This is the one change to the above code!
    start_index = 1;
    final_index = 20000;
    time = (start_index:1:final_index)/2000;
    
    scaling = 25;
    plot(time,convolution(start_index:1:final_index)*scaling);
    xt = get(gca, 'XTick');                               
    %set(gca, 'XTick', xt, 'XTickLabel', xt*1000)
    title(['contrast =',' ',num2str(contrast)]);
    width=900;
    height=600;
    set(gcf,'position',[0,0,width,height])
end


%% Part 2(b) a better attempt 



contrast_trial = 9;  %adjust this parameter manually for each trial
sigma = 20;         %tune the width

time = 1:numel(stimulus); 

smooth_fire = spike_rate_from1c(contrast_trial,sigma);  %call function (smooth firing rate from 1c)
smooth_fire = smooth_fire(1:numel(stimulus));           

trigger_curve_2a = trigger_intensity{contrast_trial};      %trigger feature intensity array for the corresponding trial
trigger_curve_2a = trigger_curve_2a(1:numel(stimulus));

trigger_curve_2a = trigger_curve_2a * max(smooth_fire)/max(trigger_curve_2a); % ~scale the trigge curve to match actual firing rate amplitude

filter_output{contrast_trial} = trigger_curve_2a;       %store trigger curve with scaling  adjustment 
adjusted_smooth_fire{contrast_trial} = smooth_fire;     %store smooth firing rate with sigma adjustment (so I don't have to write down the sigmas for 2c)
save('variables.mat','adjusted_smooth_fire','filter_output')  %save the values for future use

trigger_curve_2a(trigger_curve_2a<0)=0;             %remove values in trigger intensity <0 (replace negative values with 0)






plot(time/2000,smooth_fire,'b','LineWidth',1)
hold on
plot(time/2000,trigger_curve_2a,'r','LineWidth',1)
xlim([10,14])

title(['Trial =',' ',num2str(contrast_trial)],'FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('Firing Rate (Hz)','FontSize',16);
legend('Actual Firing Rate','Trigger Feature Intensity')
width=900;
height=600;
set(gcf,'position',[0,0,width,height])


%% Part 2(c)
for contrast_trial = 1:9
    filter = filter_output{contrast_trial};
    firing_rate = adjusted_smooth_fire{contrast_trial};
    
    subplot(3,3,contrast_trial)
    
    scatter(filter(1:10:end),firing_rate(1:10:end))
    title(['Trial =',' ',num2str(contrast_trial)],'FontSize',16);
    xlabel('Trigger Feature Intensity','FontSize',16);
    ylabel('Firing Rate (Hz)','FontSize',16);
    width=900;
    height=600;
    set(gcf,'position',[0,0,width,height])
    
end

%% Part 2(d)
for contrast_trial = 1:9
    filter = transpose(filter_output{contrast_trial});
    firing_rate = (adjusted_smooth_fire{contrast_trial});


    [N,edges] = histcounts(filter);

    for i=1:numel(edges)-1
        edge_lower = edges(i);
        edge_upper = edges(i+1);
        binned_indices = find((filter >= edge_lower & filter <=edge_upper));
        corresponding_firing_rates = firing_rate(binned_indices);
        avg = mean(corresponding_firing_rates);
        avg_firing_rate(i) = avg;

    end


    midpt = (edges(2) - edges(1))/2;
    filter_vals = edges - midpt;
    filter_vals = filter_vals(1:(numel(edges)-1));

    filter_vals = filter_vals(1:3:end);
    avg_firing_rate = avg_firing_rate(1:3:end);

    filter_vals_trial{contrast_trial} = filter_vals;
    avg_firing_rate_trial{contrast_trial} = avg_firing_rate;
    
    

    subplot(3,3,contrast_trial)
    scatter(filter_vals,avg_firing_rate)
    title(['Trial =',' ',num2str(contrast_trial)],'FontSize',16);
    xlabel('Trigger Feature Intensity','FontSize',16);
    ylabel('Firing Rate (Hz)','FontSize',16);
    width=900;
    height=600;
    set(gcf,'position',[0,0,width,height])
    
end


%% Part 2(e)
for contrast_trial = 6

    predict_fire = interp1(filter_vals_trial{contrast_trial},avg_firing_rate_trial{contrast_trial},filter_output{contrast_trial});
    
    %subplot(3,3,contrast_trial)
    plot(1:numel(predict_fire),predict_fire,'r','LineWidth',1)
    hold on
    plot(1:numel(predict_fire),adjusted_smooth_fire{contrast_trial},'b','LineWidth',1)
    xlim([10000,20000])
    xt = get(gca, 'XTick');                                
    set(gca, 'XTick', xt, 'XTickLabel', xt/2000)
    title(['Trial =',' ',num2str(contrast_trial)],'FontSize',16);
    xlabel('Time (s)','FontSize',16);
    ylabel('Firing Rate (Hz)','FontSize',16);
    legend('Predicted','Actual')
    width=900;
    height=600;
    set(gcf,'position',[0,0,width,height])

end



%% Functions

%Let's define a function that will helpn us with 2(b) which relies off of
%1c...

function out = spike_rate_from1c(contrast_trial,sigma)  %returns smooth firing rate array from 1c with length of 200,000
load('contrast_response.mat')

x = -5000:.5:5000;

gaussian = gaussmf(x,[sigma 0]);
gaussian = 2000*gaussian/(sum(gaussian));

trial = contrast_trial;
spikeTimesSelected = spikeTimes{trial};
spike_array = zeros(1,numel(stimulus));
spike_array(spikeTimesSelected)=1;
conv_firing = conv(spike_array,gaussian,'same'); %trial elemnt of array is convolution values. We will use this later!
   

[out,] = (conv_firing);

end

