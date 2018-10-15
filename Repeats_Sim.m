function Repeats_Sim

%This function runs a simulation demonstrating the effect of repeated stimuli on the frequency of illusory conjunctions.

close all

%Simulation parameters
%Temporal variables for the simulation; units are ms
Incidence_window = 5;
Sim_time = 10000;

%Define conduction delays to output neurons
Delays_1 = [6; 12]; %conduction delays to right-angle neuron
Delays_2 = [14; 22; 4]; %Conduction delays to arrow neuron

% ** Iterate through different number of distractors **

Distance = -2;
Jitter = 3; %standard deviation of the spike time

Precise_spikes = [6; 0; 12]; 
Precise_spikes(3) = Precise_spikes(3) + Distance; %Adds distance effect

Num_stim = 3; %The number of identical objects repeated that might result in an illusory conjunction

%Iterate over different numbers of objects
Percentages = zeros(3, Num_stim);
for jj = 1:Num_stim
    
    Output_count = zeros(3, 1);
    Output_num = 0;

    %Simulate each trail over the simulation time
    for ii = 1:10:Sim_time %iterate over Sim_time at 10ms intervals

        %Randomly generate Gaussean jitter to apply to each neuron
        Jitter_1 = normrnd(0, Jitter, 2, 1);
        Jitter_2 = normrnd(0, Jitter, 3, jj);

        %Calculate the arrival times of each spike, including jitter and axonal
        %conduction delay
        Arrival_times_1 = Precise_spikes(1:2) + Delays_1 + Jitter_1; %Note third diagonal neuron is excluded, as does not share connection to right angle neuron
        %Varying random jitter is added for each of the duplicated stimuli
        for kk = 1:jj
            Arrival_times_2(:, kk) = Precise_spikes + Delays_2 + Jitter_2(:, kk);
            Arrival_times_2(3, kk) = Arrival_times_2(3, kk) - 2*(kk-1); %Adds distance to the duplicate, accounting for non-uniform similarit to the first stimulus
        end
        
        %If the range of arrival times (min:max) is within the incidence
        %window, then the higher-level neuron 'activates'; if neither, count a 'non-activation'
        Incidence_1 = range(Arrival_times_1);
        Incidence_2 = range(Arrival_times_2); %returns a row vector with the range of each column

        if Incidence_1 < Incidence_window && any(Incidence_2 < Incidence_window)
            %disp('Overlap of spikes, consider increasing difference between output neurons');
            %Output_count(4) = Output_count(4) + 1;
            Output_count(1) = Output_count(1) + 0.5;
            Output_count(2) = Output_count(2) + 0.5;
        elseif Incidence_1 < Incidence_window
            Output_count(1) = Output_count(1) + 1;
        elseif any(Incidence_2 < Incidence_window)
            Output_count(2) = Output_count(2) + 1;
        else
            Output_count(3) = Output_count(3) + 1;
        end

        Output_num = Output_num + 1; 
    end


    Percentages(:, jj) = (Output_count/Output_num) * 100;
end

Percentages(2,3)/Percentages(2,1)

figure('position', [0, 0, 800, 600])
bar(1:Num_stim, Percentages(2, :))
xlabel('Number of Stimuli')
ylabel('Percentage')
