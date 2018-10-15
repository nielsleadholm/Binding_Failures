function Distance_Sim

%This function runs a simulation demonstrating the effect of spatial disparity on the frequency of illusory conjunctions.

close all

%Simulation parameters
%Temporal variables for the simulation; units are ms
Incidence_window = 5;
Sim_time = 10000;
Jitter = 3; %standard deviation of the spike time

%Define conduction delays to output neurons
Delays_1 = [6; 12]; %conduction delays to right-angle neuron
Delays_2 = [14; 22; 4]; %Conduction delays to arrow neuron
    

% ** Iterate through different amounts of spatial disparity between the objects on the
% retina **
Distance_vec = -0.5*[0:10];
Percentages = zeros(3, length(Distance_vec));
for jj = 1:length(Distance_vec)
    
    Output_count = zeros(3, 1);
    Output_num = 0;

    Precise_spikes = [6; 0; 12]; 
    Precise_spikes(3) = Precise_spikes(3) + Distance_vec(jj);

    %Simulate a cycling time series of spikes by the input neurons (+noise)
    for ii = 1:10:Sim_time %iterate over Sim_time at 10ms intervals

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

%Reverse the vectors for improved interprability
fliplr(Percentages);
Distance_vec = Distance_vec * (-1);

figure('position', [0, 0, 800, 600])
plot(Distance_vec, Percentages(1, :), 'k', 'linewidth', 1.5)
hold on
plot(Distance_vec, Percentages(2, :), 'k:', 'linewidth', 1.5)
legend('Correct Percepts', 'Illusory Conjunctions')
xlabel('Spatial Disparity')
ylabel('Percentage')
axis([0 max(Distance_vec) 0 100])
