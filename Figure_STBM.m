function Figure_STBM
% Figure for the paper "Brightness change is optimal stimulus for 
% magnocellular-projecting retinal ganglion cells" by Artem Pinchuk
% 
clc;
clear all;
load('elife-38841-fig4-data1-v2.mat')  % data from
% "Receptive field center-surround interactions mediate context-dependent 
% spatial contrast encoding in the retina"
% Maxwell H Turner, Gregory W Schwartz, Fred Rieke
% DOI: https://doi.org/10.7554/eLife.38841 

%frequency after which to cut off filter spectrum. The stimulus updated at 30
%Hz, so the cutoff should be below that value to avoid ringing
cellInd = 2;  % 1-8 Off-center   9 On-center

help_object = CenterSurroundWhiteNoise{cellInd}.stimulus.center; %center
ho_size = size(help_object);
ho_size = ho_size(2);
new_size = ho_size; % 60000; %600000; % 8000 with rate 1e4 Hz gives 0.8 sec
stimulus = zeros(new_size, 1);
for i = 1:new_size
    stimulus(i) = help_object(i);
end
help_object = CenterSurroundWhiteNoise{cellInd}.stimulus.surround; % surround
surround_stimulus = zeros(new_size, 1);
for i = 1:new_size
    surround_stimulus(i) = help_object(i);
end
% stimulus random changed each 333 ms
help_object = CenterSurroundWhiteNoise{cellInd}.response.center; %centerSurround center surround  CenterSurroundWhiteNoise{1,2}.response.centerSurround
response = zeros(new_size, 1);
for i = 1:new_size
    response(i) = help_object(i);
end
%% hi_pass_filtering
filtered_responce = hi_pass_filtering(response);
% figure(1);
% hold on
% plot(response,'r'); % responce from dataset
% plot(filtered_responce,'b');  % filtered responce 
% plot(stimulus,'black');
% hold off
%% spikes
threshold = 0.9;
spike_response = spike_detection(filtered_responce, threshold);
% figure(2);
% hold on
% plot(stimulus,'r');
% plot(spike_response,'b');
% hold off
%% STA
STA_time_period = 2000; % 2000 = 200 ms because sampleRate = 10000 Hz  
compression_coef = 100;
c_stimulus = compress_stimuly(stimulus, compression_coef);
c_spike_response = compress_spikes(spike_response, compression_coef);
ho_size = size(c_spike_response);
new_size = ho_size(1); %    1195200 
STA_c_time_period = floor(STA_time_period / compression_coef); 
STA_c = STA_calc(c_stimulus, c_spike_response, STA_c_time_period);
% figure('Name','Spike-triggered average','NumberTitle','off','Color', 'w');
% hold on
% reverse_i = zeros(STA_c_time_period, 1);
% for i = -1:-1:-STA_c_time_period
%     reverse_i(abs(i)) = i;
% end
% 
% 
% plot(reverse_i, STA_c, '--rs','LineWidth',2);
% grid on;

%grid(gca,'minor');


total_response = 0;
for frame = STA_time_period+1:new_size %
    if c_spike_response(frame) > 0.9
        total_response = total_response + 1; % count spikes
    end
end

 
alphabet = ['A','B','C','D','E','F','G','H','I','J','K' ];
empty_space = '                                                                                 ';

figure('Name', 'Spike-triggered behavior map', 'NumberTitle', 'off', 'Color', 'w', 'units','normalized','outerposition', [0 0 1 1]);

 
for time_point = 1:12-1 %1:9-1
    STA_dist_x = zeros(total_response,1);
    STA_dist_y = zeros(total_response,1);
    index = 1;
    above_diag = 0;
    on_diag = 0;
    under_diag = 0;
    for frame = STA_time_period+1:new_size 
        if c_spike_response(frame) > 0.9
            STA_dist_x(index) = c_stimulus(frame - time_point);
            STA_dist_y(index) = c_stimulus(frame - time_point - 1);
            if STA_dist_y(index) > STA_dist_x(index)
                above_diag = above_diag + 1;
            end
            if STA_dist_y(index) == STA_dist_x(index)
                on_diag = on_diag + 1;
            end
            if STA_dist_y(index) < STA_dist_x(index)
                under_diag = under_diag + 1;
            end
            
            index = index +1; 
        end
    end
	STA_x = STA_c(time_point);
	STA_y = STA_c(time_point + 1);
%    subplot(3,3,time_point);
    subplot(4,3,time_point);
    hold on
    plot(STA_dist_x,STA_dist_y,'.','MarkerFaceColor','g','MarkerSize',10);
    title({['\bf\fontsize{11}', alphabet(time_point), empty_space];['\rm\fontsize{10}Points above/on/under diagonal = ',num2str(above_diag), '/',num2str(on_diag), '/',num2str(under_diag)]});
    %title(['Points above/on/under diagonal = ',num2str(above_diag), '/',num2str(on_diag), '/',num2str(under_diag)]);
    plot(STA_x,STA_y,'rs','MarkerEdgeColor','red','MarkerSize',9,'LineWidth',3);
    plot([-1.5,1.5],[0,0],'-red');
    plot([0,0],[-1.5,1.5],'-black');
    %axis([-1 1 -1 1 ]);
    
    set(gca,'XTick',-0.5:0.5:0.5,'YTick',-0.5:0.5:0.5);
    axis([-0.75 0.75 -0.75 0.75 ]);
    xlabel(['Stimuli for t = ' num2str(-time_point*10) ' ms']);
    ylabel(['Stimuli for t = ' num2str(-(time_point + 1)*10) ' ms']);
    hold off
end

subplot(4,3,12, 'YAxisLocation', 'right');
hold on
    reverse_i = zeros(STA_c_time_period, 1);
    for i = -1:-1:-STA_c_time_period
        reverse_i(abs(i)) = 10*i;
    end

    title({['\bf\fontsize{11}L', empty_space];'\rm\fontsize{10}STA'}); %,'HorizontalAlignment','left'

    %title({'First line';'Second line'},'HorizontalAlignment','left')
    xlabel('Time before spike [ms]');
    plot(reverse_i, STA_c, '-rs','LineWidth',2 );
    axis([-200 0 -0.4 0.2 ]);
    % axes('YAxisLocation', 'left');

    %grid(gca,'minor');
    grid on;
hold off




% Defaults for this blog post
width = 4*4.4;     % Width in inches for subplot(3,4,time_point);
height = 4*3;    % Height in inches
width = 4*3;     % Width in inches for subplot(4,3,time_point);
height = 4*4;    % Height in inches

set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'InvertHardcopy','off');
set(gcf,'PaperUnits', 'inches');
%papersize = get(gcf, 'PaperSize');
left = 0;%(papersize(1)- width)/2;
bottom = 0;%(papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

% Save the file as PNG
print('Figure_STBM  LOW_RES','-dpng','-r100');
print('Figure_STBM','-dpng','-r300');

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