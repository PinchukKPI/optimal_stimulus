function Figure_data_observation()
% Figure for the paper "Brightness change is optimal stimulus for 
% magnocellular-projecting retinal ganglion cells" by Artem Pinchuk
% 
% artem.pinchuk at gmail.com
% 20.10.2021
clc;
clear all;

load('elife-38841-fig4-data1-v2.mat')  % source data from
% "Receptive field center-surround interactions mediate context-dependent 
% spatial contrast encoding in the retina"
% Maxwell H Turner, Gregory W Schwartz, Fred Rieke
% DOI: https://doi.org/10.7554/eLife.38841 
% https://elifesciences.org/articles/38841
% CenterSurroundWhiteNoise is a cell array, each entry of which is a structure that corresponds
% to a cell in the dataset. Structure fields:
% .stimulus = concatenated stimulus traces for .center and .surround
% stimuli
% .response = concatenated excitatory conductance response traces (in nS)
% for .center, .surround and .centerSurround stimuli
% Note that data are concatenated and grouped for convenience, but were
% acquired in interleaved trials


 
OFF_cellInd = 2;  % 1-8 Off-cell
ON_cellInd = 13;  %   9-15 On-cell
STA_time_period = 2000; 
%% read OFF cell data
help_object = CenterSurroundWhiteNoise{OFF_cellInd}.stimulus.center; %center
ho_size = size(help_object);
ho_size = ho_size(2);
new_size = ho_size; % 60000; %600000; % 8000 with rate 1e4 Hz gives 0.8 sec

OFF_center_stimulus = zeros(new_size, 1);
for i = 1:new_size
    OFF_center_stimulus(i) = help_object(i);% stimulus random changed each 33 ms?
end
help_object = CenterSurroundWhiteNoise{OFF_cellInd}.stimulus.surround; % surround
OFF_surround_stimulus = zeros(new_size, 1);
for i = 1:new_size
    OFF_surround_stimulus(i) = help_object(i);
end

help_object = CenterSurroundWhiteNoise{OFF_cellInd}.response.center; %centerSurround center surround  
OFF_center_response = zeros(new_size, 1);
for i = 1:new_size
    OFF_center_response(i) = help_object(i);
end
help_object = CenterSurroundWhiteNoise{OFF_cellInd}.response.surround; %centerSurround center surround  
OFF_surround_response = zeros(new_size, 1);
for i = 1:new_size
    OFF_surround_response(i) = help_object(i);
end
%% read ON cell data

help_object = CenterSurroundWhiteNoise{ON_cellInd}.stimulus.center;
ho_size = size(help_object);
ho_size = ho_size(2);
new_size = ho_size;
ON_center_stimulus = zeros(new_size, 1);
for i = 1:new_size
    ON_center_stimulus(i) = help_object(i);
end
help_object = CenterSurroundWhiteNoise{ON_cellInd}.stimulus.surround; 
ON_surround_stimulus = zeros(new_size, 1);
for i = 1:new_size
    ON_surround_stimulus(i) = help_object(i);
end
% stimulus random changed each 333 ms
help_object = CenterSurroundWhiteNoise{ON_cellInd}.response.center; 
ON_center_response = zeros(new_size, 1);
for i = 1:new_size
    ON_center_response(i) = help_object(i);
end
help_object = CenterSurroundWhiteNoise{ON_cellInd}.response.surround;
ON_surround_response = zeros(new_size, 1);
for i = 1:new_size
    ON_surround_response(i) = help_object(i);
end

%% hi_pass_filtering
OFF_center_filtered_responce = hi_pass_filtering(OFF_center_response);
OFF_surround_filtered_responce = hi_pass_filtering(OFF_surround_response);
ON_center_filtered_responce = hi_pass_filtering(ON_center_response);
ON_surround_filtered_responce = hi_pass_filtering(ON_surround_response);

% figure(1);
% hold on
% plot(OFF_center_response,'r'); % responce from dataset
% plot(OFF_center_stimulus,'black');
% hold off
%% spikes
threshold = 0.9;
OFF_center_spike_response = spike_detection(OFF_center_filtered_responce, threshold);
OFF_surround_spike_response = spike_detection(OFF_surround_filtered_responce, threshold);
ON_center_spike_response = spike_detection(ON_center_filtered_responce, threshold);
ON_surround_spike_response = spike_detection(ON_surround_filtered_responce, threshold);

% figure(2);
% hold on
% plot(OFF_center_filtered_responce,'r');
% plot(OFF_center_spike_response,'b');
% hold off

%% OFF CENTRAL find spikes for removing opposite stimulus for central region
OFF_center_spike_frames = find_spike_frames(OFF_center_spike_response, STA_time_period);
OFF_center_spikes_number = find_spikes_number(OFF_center_spike_response, STA_time_period);
FLAG_surround = 1;
FLAG_center = 0;
disp('Spike in case of removing opposite stimulus for "center region" of parasol cell');
find_and_display_opposite_stimuli(FLAG_center, OFF_cellInd, OFF_center_spikes_number, OFF_center_spike_frames, OFF_center_stimulus, STA_time_period);

%% OFF SURROUND find spikes for removing opposite stimulus for surround region
OFF_surround_spike_frames = find_spike_frames(OFF_surround_spike_response, STA_time_period);
OFF_surround_spikes_number = find_spikes_number(OFF_surround_spike_response, STA_time_period);
disp(' ');
disp('Spike in case of removing opposite stimulus for "surround region" of parasol cell');
find_and_display_opposite_stimuli(FLAG_surround, OFF_cellInd, OFF_surround_spikes_number, OFF_surround_spike_frames, OFF_surround_stimulus, STA_time_period);



%% ON CENTRAL find spikes for removing opposite stimulus for central region
ON_center_spike_frames = find_spike_frames(ON_center_spike_response, STA_time_period);
ON_center_spikes_number = find_spikes_number(ON_center_spike_response, STA_time_period);
disp(' ');
disp('Spike in case of removing opposite stimulus for "center region" of parasol cell');
find_and_display_opposite_stimuli(FLAG_center, ON_cellInd, ON_center_spikes_number, ON_center_spike_frames, ON_center_stimulus, STA_time_period);

%% ON SURROUND find spikes for removing opposite stimulus for surround region
ON_surround_spike_frames = find_spike_frames(ON_surround_spike_response, STA_time_period);
ON_surround_spikes_number = find_spikes_number(ON_surround_spike_response, STA_time_period);
disp(' ');
disp('Spike in case of removing opposite stimulus for "surround region" of parasol cell');
find_and_display_opposite_stimuli(FLAG_surround, ON_cellInd, ON_surround_spikes_number, ON_surround_spike_frames, ON_surround_stimulus, STA_time_period);

% we choose to show following sequence of stimuli which trigger spikes 
ON_S_show_spike_stimulus = 423;
ON_C_show_spike_stimulus = 1333;
OFF_S_show_spike_stimulus = 475;
OFF_C_show_spike_stimulus = 3418;

%% count spikes
 
STA_OFF_C = STA_calc(OFF_center_stimulus, OFF_center_spike_response, STA_time_period);
STA_OFF_S = STA_calc(OFF_surround_stimulus, OFF_surround_spike_response, STA_time_period);
STA_ON_C = STA_calc(ON_center_stimulus, ON_center_spike_response, STA_time_period);
STA_ON_S = STA_calc(ON_surround_stimulus, ON_surround_spike_response, STA_time_period);

OFF_C_all_spike_triggered = zeros(STA_time_period,OFF_center_spikes_number);
for spike = 1:OFF_center_spikes_number
    frame = OFF_center_spike_frames(spike);
    for i = 1:STA_time_period
        OFF_C_all_spike_triggered(i,spike) =  OFF_center_stimulus(frame - i);
    end
end

OFF_S_all_spike_triggered = zeros(STA_time_period,OFF_surround_spikes_number);
for spike = 1:OFF_surround_spikes_number
    frame = OFF_surround_spike_frames(spike); 
    for i = 1:STA_time_period
        OFF_S_all_spike_triggered(i,spike) =  OFF_surround_stimulus(frame - i);
    end
end

ON_C_all_spike_triggered = zeros(STA_time_period,ON_center_spikes_number);
for spike = 1:ON_center_spikes_number
    frame = ON_center_spike_frames(spike);
    for i = 1:STA_time_period
        ON_C_all_spike_triggered(i,spike) =  ON_center_stimulus(frame - i);
    end
end

ON_S_all_spike_triggered = zeros(STA_time_period,ON_surround_spikes_number);
for spike = 1:ON_surround_spikes_number
    frame = ON_surround_spike_frames(spike); 
    for i = 1:STA_time_period
        ON_S_all_spike_triggered(i,spike) =  ON_surround_stimulus(frame - i);
    end
end

 
alphabet = ['A','B','C','D','E','F','G','H','I','J','K' ];
empty_space   = '                                                                               ';

figure('Name', 'Data observation','NumberTitle','off','Color', 'w', 'units','inches'); 
clf;
%set(gcf,'Resize','off');



subplot(2,2,1, 'YAxisLocation', 'right'); hold on;
hold on
plot( -(1:STA_time_period)*0.1, OFF_C_all_spike_triggered(:,OFF_C_show_spike_stimulus), '-b','LineWidth',1.5 );
plot( -(1:STA_time_period)*0.1, STA_OFF_C,'r', 'LineWidth',1.5);
plot( -(1:STA_time_period)*0.1, zeros(STA_time_period, 1),'Color','black','LineWidth',1);
%grid on
axis([-STA_time_period*0.1 0 -1 1 ]);
	title({['\bf\fontsize{11}', alphabet(1), empty_space];['\rm\fontsize{10} OFF parasol cell']});
    legend('Opposite stimulus','Center STA', 'Location','NorthWest');
    xlabel('Time [ms]');
    ylabel('Brightness');
    %set(gca,'YTickLabel',{' '});
    lh = legend();
    set(lh, 'Box', 'off');
    set(lh, 'Color', 'none');
hold off


subplot(2,2,2, 'YAxisLocation', 'right'); hold on;
hold on
plot( -(1:STA_time_period)*0.1, OFF_S_all_spike_triggered(:,OFF_S_show_spike_stimulus), '-b','LineWidth',1.5 );   % 475
plot( -(1:STA_time_period)*0.1, STA_OFF_S,'r', 'LineWidth',1.5);
plot( -(1:STA_time_period)*0.1, zeros(STA_time_period, 1),'Color','black','LineWidth',1);
%grid on
axis([-STA_time_period*0.1 0 -1 1 ]);
	title({['\bf\fontsize{11}', alphabet(2), empty_space];['\rm\fontsize{10} OFF parasol cell']});
    legend('Opposite stimulus','Surround STA', 'Location','NorthWest');
    xlabel('Time [ms]');
    ylabel('Brightness');
    %set(gca,'YTickLabel',{' '});
    lh = legend();
    set(lh, 'Box', 'off');
    set(lh, 'Color', 'none');
hold off

 %% ON Cell plot
subplot(2,2,3, 'YAxisLocation', 'right'); hold on;
hold on
plot( -(1:STA_time_period)*0.1, ON_C_all_spike_triggered(:,ON_C_show_spike_stimulus), '-b','LineWidth',1.5 );
plot( -(1:STA_time_period)*0.1, 1.3*STA_ON_C,'r', 'LineWidth',1.5);
plot( -(1:STA_time_period)*0.1, zeros(STA_time_period, 1),'Color','black','LineWidth',1);
%grid on
axis([-STA_time_period*0.1 0 -1 1 ]);
	title({['\bf\fontsize{11}', alphabet(3), empty_space];['\rm\fontsize{10} ON parasol cell']});
    legend('Opposite stimulus','Center STA', 'Location','NorthWest');
    xlabel('Time [ms]');
    ylabel('Brightness');
    %set(gca,'YTickLabel',{' '});
    lh = legend();
    set(lh, 'Box', 'off');
    set(lh, 'Color', 'none');
hold off


subplot(2,2,4, 'YAxisLocation', 'right'); hold on;
hold on
plot( -(1:STA_time_period)*0.1, ON_S_all_spike_triggered(:,ON_S_show_spike_stimulus), '-b','LineWidth',1.5 );   % 423 for #13
plot( -(1:STA_time_period)*0.1, 1.3*STA_ON_S,'r', 'LineWidth',1.5);
plot( -(1:STA_time_period)*0.1, zeros(STA_time_period, 1),'Color','black','LineWidth',1);
%grid on
axis([-STA_time_period*0.1 0 -1 1 ]);
	title({['\bf\fontsize{11}', alphabet(4), empty_space];['\rm\fontsize{10} ON parasol cell']});
    legend('Opposite stimulus','Surround STA', 'Location','NorthWest');
    xlabel('Time [ms]');
    ylabel('Brightness');
    %set(gca,'YTickLabel',{' '});
    lh = legend();
    set(lh, 'Box', 'off');
    set(lh, 'Color', 'none');
hold off

% Defaults for this blog post
width = 3*3;     % Width in inches 
height = 3*3;    % Height in inches

set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
% Here we preserve the size of the image when we save it.
%set(gcf,'InvertHardcopy','on');
set(gcf,'InvertHardcopy','off');
set(gcf,'PaperUnits', 'inches');
%papersize = get(gcf, 'PaperSize');
left = 0;%(papersize(1)- width)/2;
bottom = 0;%(papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
 
set(gcf,'Position', myfiguresize);
set(gcf,'Resize','off');
% Save the file as PNG
print('Figure_data_observation','-dpng','-r300');
print('Figure_data_observation LOW_RES','-dpng','-r100');


end

%%  Support functions
function filtered = hi_pass_filtering(input_array)
    new_size = size(input_array);
    new_size = new_size(1);

     % 5 '-0.1'  1 '1'  5 '-0.1'  = 11
    kernel = [ -0.1 -0.1 -0.1 -0.1 -0.1  1  -0.1 -0.1 -0.1 -0.1 -0.1 ];
    kernel_size = 11;

    filtered = zeros(new_size, 1) ;
    for i = 1:new_size
        kernel_sum = 0;
        for filter_index = 1:kernel_size
            if (i + filter_index - 6) > 0 && (i + filter_index - 6) < new_size
                kernel_sum = kernel_sum + kernel(filter_index)* input_array(i + filter_index - 6);
            end
            
        end
        filtered(i) = kernel_sum;
    end
end


function find_and_display_opposite_stimuli(FLAG, cellInd, spikes_number, spike_frames, stimulus, STA_time_period)
% FLAG = 1; FLAG_surround
% FLAG = 0; FLAG_center

% 1:332 one stimulus duration window  1-332  333-664 665-996  (332)
% for spike triggered sequence of stimuli inside 332-664 664-996 we find 2 changes of stimulus value

% all 3 steps > 0  and step(1) < step(2) < step(3)  for OFF-cell center
% region and for On-cell surround region
    for spike = 1:spikes_number
        frame = spike_frames(spike);
        if frame > STA_time_period
            for i = 664:996
                if stimulus(frame - i) < stimulus(frame - i - 1) && stimulus(frame - i) > 0
                    for e = 332:664
                        if stimulus(frame - e) < stimulus(frame - e - 1) && stimulus(frame - e) > 0
                            if FLAG == 0 && cellInd < 9 % 1-8 Off-cell center   9 On-center
                                disp(['OFF-cell #',num2str(cellInd),' "center region" spike#=',num2str(spike), ' frame start from ',num2str(frame) ]);
                            end
                            if FLAG == 1 && cellInd >= 9  % On-cell surround
                                disp(['ON-cell #',num2str(cellInd),' "surround region" spike#=',num2str(spike), ' frame start from ',num2str(frame) ]);
                            end
                        end
                    end
                end
                % all 3 steps < 0  and step(1) > step(2) > step(3)  for On-cell center
                % region and Off-cell surround region
                if stimulus(frame - i) > stimulus(frame - i - 1) && stimulus(frame - i) < 0
                    for e = 332:664
                        if stimulus(frame - e) > stimulus(frame - e - 1) && stimulus(frame - e) < 0
                            if FLAG == 0 && cellInd >= 9 % On-cell center
                                disp(['ON-cell #',num2str(cellInd),' "center region" spike#=',num2str(spike), ' frame start from ',num2str(frame) ]);
                            end
                            if FLAG == 1 && cellInd < 9  % Off-cell surround
                                disp(['OFF-cell #',num2str(cellInd),' "surround region" spike#=',num2str(spike), ' frame start from ',num2str(frame) ]);
                            end
                        end
                    end
                end
            end
        end
    end
end

function spike_frames = find_spike_frames(spike_response, STA_time_period)
    spikes_length = length(spike_response);
    spikes_number = 0;
    for i = STA_time_period:spikes_length
        if spike_response(i) == 1
            spikes_number = spikes_number + 1;
        end
    end

    spike_frames = zeros(spikes_number,1);
    spikes_number = 0;
    for i = STA_time_period:spikes_length
        if spike_response(i) == 1
            spikes_number = spikes_number + 1;
            spike_frames(spikes_number) = i;
        end
    end
end

function spikes_number = find_spikes_number(spike_response, STA_time_period)
    spikes_length = length(spike_response);
    spikes_number = 0;
    for i = STA_time_period:spikes_length
        if spike_response(i) == 1
            spikes_number = spikes_number + 1;
        end
    end
end

function spikes = spike_detection(input_array, threshold)
    recived_size = size(input_array);
    spikes = zeros(recived_size(1), 1);
    for i = 1:recived_size - 1
        if input_array(i + 1) <= threshold && input_array(i) > threshold
            spikes(i) = 1;
        end
    end
end

function STA = STA_calc(input_stimulus, input_spike_response, STA_time_period )
    total_response = 0;
    new_size = size(input_spike_response);
    STA = zeros(STA_time_period,1);
    for frame = STA_time_period+1:new_size %  
        if input_spike_response(frame) > 0.9
            total_response = total_response + 1;
            for i = 1:STA_time_period
                STA(i) =  STA(i) +  input_stimulus(frame - i);
            end
        end
    end
    STA  =  STA / total_response;
end

function compressed = compress_stimuly(input_stimulus, compression_coef)
% return just last point of each compressed part of stimulus
% do not use for stimuly difference
    recived_size = size(input_stimulus);
    compressed_size = floor(recived_size(1) / compression_coef);
    compressed = zeros(compressed_size, 1);
    for i = 1:compressed_size
        compressed(i) = input_stimulus(i*compression_coef);
    end
end

function compressed = compress_spikes(input_spike, compression_coef)
% can be used to compress stiumulus difference until stimulus rate size
    recived_size = size(input_spike);
    compressed_size = floor(recived_size(1) / compression_coef);
    compressed = zeros(compressed_size, 1);
    for i = 1:compressed_size
        spikes_counted = 0;
        for a = 1:compression_coef
            if ((i-1)*compression_coef + a) <= recived_size(1)
                spikes_counted = spikes_counted + input_spike((i-1)*compression_coef + a);
            end
        end
        compressed(i) = spikes_counted;
    end
end


% Spike in case of removing opposite stimulus for "center region" of parasol cell
% OFF-cell #2 "center region" spike#=394 frame start from 68652
% OFF-cell #2 "center region" spike#=395 frame start from 68695
% OFF-cell #2 "center region" spike#=445 frame start from 77745
% OFF-cell #2 "center region" spike#=446 frame start from 77822
% OFF-cell #2 "center region" spike#=447 frame start from 77915
% OFF-cell #2 "center region" spike#=448 frame start from 77928
% OFF-cell #2 "center region" spike#=449 frame start from 77971
% OFF-cell #2 "center region" spike#=450 frame start from 78123
% OFF-cell #2 "center region" spike#=451 frame start from 78154
% OFF-cell #2 "center region" spike#=799 frame start from 155373
% OFF-cell #2 "center region" spike#=1272 frame start from 265335
% OFF-cell #2 "center region" spike#=1677 frame start from 372724
% OFF-cell #2 "center region" spike#=1678 frame start from 372788
% OFF-cell #2 "center region" spike#=1679 frame start from 372807
% OFF-cell #2 "center region" spike#=1718 frame start from 383461
% OFF-cell #2 "center region" spike#=2388 frame start from 522031
% OFF-cell #2 "center region" spike#=2422 frame start from 532627
% OFF-cell #2 "center region" spike#=2560 frame start from 556213
% OFF-cell #2 "center region" spike#=2561 frame start from 556343
% OFF-cell #2 "center region" spike#=2562 frame start from 556368
% OFF-cell #2 "center region" spike#=2563 frame start from 556380
% OFF-cell #2 "center region" spike#=2564 frame start from 556417
% OFF-cell #2 "center region" spike#=2659 frame start from 570561
% OFF-cell #2 "center region" spike#=2660 frame start from 570677
% OFF-cell #2 "center region" spike#=2721 frame start from 587382
% OFF-cell #2 "center region" spike#=2722 frame start from 587401
% OFF-cell #2 "center region" spike#=2723 frame start from 587514
% OFF-cell #2 "center region" spike#=2724 frame start from 587526
% OFF-cell #2 "center region" spike#=2899 frame start from 631797
% OFF-cell #2 "center region" spike#=3054 frame start from 662553
% OFF-cell #2 "center region" spike#=3055 frame start from 662583
% OFF-cell #2 "center region" spike#=3067 frame start from 664902
% OFF-cell #2 "center region" spike#=3418 frame start from 739596  perfect
% OFF-cell #2 "center region" spike#=3751 frame start from 1088297
% not lower than 0


% OFF-cell #2 "surround region" spike#=16 frame start from 8564
% OFF-cell #2 "surround region" spike#=17 frame start from 8579
% OFF-cell #2 "surround region" spike#=18 frame start from 8613
% OFF-cell #2 "surround region" spike#=172 frame start from 58399
% OFF-cell #2 "surround region" spike#=474 frame start from 144695
% OFF-cell #2 "surround region" spike#=475 frame start from 144721 perfect
% OFF-cell #2 "surround region" spike#=841 frame start from 281799
% OFF-cell #2 "surround region" spike#=842 frame start from 281846
% OFF-cell #2 "surround region" spike#=1786 frame start from 572310
% OFF-cell #2 "surround region" spike#=1787 frame start from 572332
% OFF-cell #2 "surround region" spike#=2097 frame start from 660391
% OFF-cell #2 "surround region" spike#=2139 frame start from 672802
% OFF-cell #2 "surround region" spike#=2140 frame start from 672812
% OFF-cell #2 "surround region" spike#=2141 frame start from 672825
% OFF-cell #2 "surround region" spike#=2474 frame start from 781151
% OFF-cell #2 "surround region" spike#=2492 frame start from 790462
% OFF-cell #2 "surround region" spike#=2536 frame start from 1025452
% 
% Spike in case of removing opposite stimulus for "center region" of parasol cell
% ON-cell #13 "center region" spike#=20 frame start from 10014
% ON-cell #13 "center region" spike#=126 frame start from 38993
% ON-cell #13 "center region" spike#=127 frame start from 39019
% ON-cell #13 "center region" spike#=261 frame start from 70679
% ON-cell #13 "center region" spike#=273 frame start from 76993
% ON-cell #13 "center region" spike#=274 frame start from 77007
% ON-cell #13 "center region" spike#=546 frame start from 146414
% ON-cell #13 "center region" spike#=559 frame start from 150521
% ON-cell #13 "center region" spike#=560 frame start from 150526
% ON-cell #13 "center region" spike#=561 frame start from 150727
% ON-cell #13 "center region" spike#=709 frame start from 191331
% ON-cell #13 "center region" spike#=848 frame start from 224101
% ON-cell #13 "center region" spike#=893 frame start from 239041
% ON-cell #13 "center region" spike#=942 frame start from 252516
% ON-cell #13 "center region" spike#=943 frame start from 252541
% ON-cell #13 "center region" spike#=944 frame start from 252771
% ON-cell #13 "center region" spike#=945 frame start from 252829
% ON-cell #13 "center region" spike#=946 frame start from 252927
% ON-cell #13 "center region" spike#=998 frame start from 269028
% ON-cell #13 "center region" spike#=999 frame start from 269180
% ON-cell #13 "center region" spike#=1000 frame start from 269245
% ON-cell #13 "center region" spike#=1051 frame start from 285843
% ON-cell #13 "center region" spike#=1052 frame start from 285939
% ON-cell #13 "center region" spike#=1053 frame start from 286162
% ON-cell #13 "center region" spike#=1086 frame start from 293450
% ON-cell #13 "center region" spike#=1137 frame start from 310103
% ON-cell #13 "center region" spike#=1138 frame start from 310122
% ON-cell #13 "center region" spike#=1139 frame start from 310261
% ON-cell #13 "center region" spike#=1140 frame start from 310269
% ON-cell #13 "center region" spike#=1141 frame start from 310293
% ON-cell #13 "center region" spike#=1142 frame start from 310346
% ON-cell #13 "center region" spike#=1186 frame start from 320963
% ON-cell #13 "center region" spike#=1249 frame start from 335186
% ON-cell #13 "center region" spike#=1286 frame start from 349374
% ON-cell #13 "center region" spike#=1287 frame start from 349522
% ON-cell #13 "center region" spike#=1297 frame start from 353534
% ON-cell #13 "center region" spike#=1298 frame start from 353543
% ON-cell #13 "center region" spike#=1332 frame start from 362009
% ON-cell #13 "center region" spike#=1333 frame start from 362034  perfect
%  
% Spike in case of removing opposite stimulus for "surround region" of parasol cell
% ON-cell #13 "surround region" spike#=423 frame start from 93856   perfect
% ON-cell #13 "surround region" spike#=424 frame start from 93922
% ON-cell #13 "surround region" spike#=558 frame start from 131490
% ON-cell #13 "surround region" spike#=559 frame start from 131564
% ON-cell #13 "surround region" spike#=560 frame start from 131585
% ON-cell #13 "surround region" spike#=561 frame start from 131782
% ON-cell #13 "surround region" spike#=892 frame start from 215241