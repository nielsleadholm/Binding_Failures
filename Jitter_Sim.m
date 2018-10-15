function Jitter_Sim

%This function runs a simulation demonstrating the effect of jitter on the frequency of illusory conjunctions.

close all

%Simulation parameters
%Temporal variables for the simulation; units are ms
Incidence_window = 5;
Sim_time = 10000;

%Input neuron's spike time
    %Base case; the below must be uncommented to use the stimuli
    %Precise_spikes = [6; 0; 2]; %Input stimuli are 'distant' on
    %the retina, and thus the 3rd neuron's time is inconsistent with an
    %illusory conjunction
    
    %Case where stimuli are very close to one-another in the periphery
    Precise_spikes = [6; 0; 12];
    

%Define conduction delays to output neurons
Delays_1 = [6; 12]; %conduction delays to right-angle neuron
Delays_2 = [14; 22; 4]; %Conduction delays to arrow neuron

% **Iterate over different values of jitter **
Jitter_vec = 0.5*[0:20];
Percentages = zeros(3, length(Jitter_vec)); %Initialise array for holding percentage frequencies
for jj = 1:length(Jitter_vec)
    Jitter = Jitter_vec(jj);
    
    Output_count = zeros(3, 1); %Tallies the number of incidences of each output possibility
    %The first and second row represent the right angle and arrow neurons
    %respectively
    %The 3rd represents no activation
    Output_num = 0; %Count total number of activations performed

    %Iterate through each trial over the simulation time
    for ii = 1:10:Sim_time 

        %Randomly generate zero-mean Gaussean jitter to apply to each neuron
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

    Percentages(:, jj) = (Output_count/Output_num) * 100;
end

Relative_percent_j = Percentages(2,:) ./ Percentages(1, :) * 100; %Relative incidence of ICs vs. correct percepts

figure('position', [0, 0, 1000, 800])
plot(Jitter_vec, Percentages(1, :), 'k', 'linewidth', 1.5)
hold on
plot(Jitter_vec, Percentages(2, :), 'k:', 'linewidth', 1.5)
plot(Jitter_vec, Percentages(3, :), 'k-o', 'linewidth', 1.5)
plot(Jitter_vec, Relative_percent_j, 'k--', 'linewidth', 1.5)
legend('Correct Percepts', 'Illusory Conjunctions', 'Inattentional Blindness', 'Relative Percent ICs')
xlabel('Jitter (ms)')
ylabel('Percentage')
axis([0 max(Jitter_vec) 0 100])
