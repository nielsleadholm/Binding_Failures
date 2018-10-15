function Fixation_Sim

close all

%This function runs a simulation demonstrating the effect of retinal
%fixation on coincidence detection in the model.

%Simulation parameters
%Temporal variables for the simulation; units are ms
Incidence_window = 5;
Sim_time = 10000;
Retinal_fixation = 1; %Retinal fixation by default is not used

Output_count = zeros(3, 1); %Tallies the number of incidences of each output possibility
%The first and second row represent the right angle and arrow neurons
%The 3rd represents no activation
Output_num = 0; %Counter for the number of trials run

%Define conduction delays to output neurons
Delays_1 = [6; 12]; %conduction delays to right-angle neuron
Delays_2 = [14; 22; 4]; %Conduction delays to arrow neuron

%Iterate through each trial over the simulation time
for ii = 1:10:Sim_time %iterate over Sim_time at 10ms intervals

    %Perform the below if retinal fixation simulation
    if Retinal_fixation == true
        Latency_max = 30;
        Precise_spikes = randi(Latency_max, 3, 1); %Creates a vector of random integer values between 1 and Latency_max
        Jitter = 0; %Removes jitter as not concerned about absence of attention
    end

    %Randomly generate Gaussean jitter to apply to each neuron
    Jitter_1 = normrnd(0, Jitter, 2, 1);
    Jitter_2 = normrnd(0, Jitter, 3, 1);

    %Calculate the arrival times of each spike, including jitter and axonal
    %conduction delay
    Arrival_times_1 = Precise_spikes(1:2) + Delays_1 + Jitter_1; %Note third diagonal neuron is excluded, as does not share connection to right angle neuron
    Arrival_times_2 = Precise_spikes + Delays_2 + Jitter_2;

    %If the range of arrival times (min:max) is within the incidence
    %window, then the higher-level neuron 'activates'; if neither, count a 'non-activation'
    Incidence_1 = range(Arrival_times_1);
    Incidence_2 = range(Arrival_times_2);

    if Incidence_1 < Incidence_window && Incidence_2 < Incidence_window
        %If signal arrives in both neurons temporal windows, then each
        %neuron registers half an activation
        Output_count(1) = Output_count(1) + 0.5;
        Output_count(2) = Output_count(2) + 0.5;
    elseif Incidence_1 < Incidence_window
        Output_count(1) = Output_count(1) + 1;
    elseif Incidence_2 < Incidence_window
        Output_count(2) = Output_count(2) + 1;
    else
        Output_count(3) = Output_count(3) + 1;
    end

    Output_num = Output_num + 1; 
end

Percentages = (Output_count/Output_num) * 100;
Bar_percents = [Percentages(1)+Percentages(2); Percentages(3)];

figure('position', [0, 0, 800, 600])
bar(Bar_percents, 'k');
l{1}='Coincident Activation'; l{2}='Non-Coincidence';    
set(gca,'xticklabel', l);
ylabel('Percentage')
