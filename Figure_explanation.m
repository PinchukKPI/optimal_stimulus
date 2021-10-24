function Figure_explanation()
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


cellInd = 14;     % 1-8 Off-center   9-15 On-center
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

help_object = CenterSurroundWhiteNoise{cellInd}.response.center; %centerSurround center surround 
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

compression_coef = 166; % stimulus changed each 332;
c_stimulus = compress_stimuly(stimulus, compression_coef);
c_spike_response = compress_spikes(spike_response, compression_coef);

%% Figure
alphabet = ['A','B','C','D','E','F','G','H','I','J','K' ];
empty_space   = '                                                                               ';
gray_color = [.45 0.45 0.45];
%figure('Name','STA as filtet','NumberTitle','off','Color', 'w', 'units','normalized','outerposition', [0 0 1 1]); 
figure('Name','Explanation','NumberTitle','off','Color', 'w', 'units','inches'); 
clf;
%set(gcf,'Resize','off');

%% stimuli distribution histogram
distribution_borders = [-0.95 -0.85 -0.75 -0.65 -0.55 -0.45 -0.35 -0.25 -0.15 -0.05 0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95];
distribution = zeros(20,1);
distr_x = zeros(20,1);
for e = 1:20
    distr_x(e) = (e - 10) * 0.1;
end
for i = 1:floor(new_size/compression_coef)
    for e = 1:19
        if c_stimulus(i) > distribution_borders(e) && c_stimulus(i) < distribution_borders(e + 1)
            distribution(e) = distribution(e) + 1;
        end
    end
end

subplot(1,3,1);

barh( distr_x, distribution,'hist'); %, 'black','LineWidth',1.5); % we have gaussian distribution 
title({['\bf\fontsize{11}', alphabet(1), empty_space]});
axis([0 600  -0.75 0.75 ]);
set(gca,'YGrid','on');
ylabel('Stimulus');
%% Stimumlus
subplot(1,3,2);
title({['\bf\fontsize{11}', alphabet(2), empty_space]});
hold on
%stem(c_spike_response(3500:4000) * 0.5,'black');
start_t  = 360;
finish_t = 780;
time_line = zeros(finish_t - start_t,1);
for i = 1:finish_t - start_t + 1
    time_line(i) = fix((i+1)/2); 
end
plot(time_line, c_stimulus(start_t:finish_t),'b', 'LineWidth',1.5);
plot(time_line(29:30), c_stimulus(start_t + 29 - 1:start_t + 30 - 1),'r','LineWidth',2.5); % optimal stimulus
plot(time_line(31:32), c_stimulus(start_t + 31 - 1:start_t + 32 - 1),'black','LineWidth',2.5); % opposite stimulus

plot(time_line(30), c_stimulus(start_t + 30 - 1), 'ro','MarkerEdgeColor','red','MarkerSize',7,'LineWidth',2.5); % optimal stimulus
plot(time_line(32), c_stimulus(start_t + 32 - 1),'ro','MarkerEdgeColor','black','MarkerSize',7,'LineWidth',2.5); % opposite stimulus
%plot(surround_stimulus(5000:10000),'r');
hold off
set(gca,'YGrid','on');
set(gca,'YLim', [-0.75, 0.75]);
set(gca,'YTick',-1:0.2:1);
%set(gca,'YTickLabel',{'0', '0' });
 


set(gca,'XLim', [10, 40]);
%set(gca,'XTick',0:250:250);
set(gca,'XTickLabel','  ');
ylabel('Stimulus');
xlabel('Time');
%axis([0 200 -0.75 0.75 ]);

%% Spike triggered behaviour map example

    subplot(1,3,3);
    hold on
    %plot(STA_dist_x,STA_dist_y,'.','MarkerFaceColor','g','MarkerSize',10);
    title({['\bf\fontsize{11}', alphabet(3), empty_space]});
    %title(['Points above/on/under diagonal = ',num2str(above_diag), '/',num2str(on_diag), '/',num2str(under_diag)]);
    %plot(STA_x,STA_y,'rs','MarkerEdgeColor','red','MarkerSize',9,'LineWidth',3);
    plot([-1.5,1.5],[0,0],'-red');
    plot([0,0],[-1.5,1.5],'-black');
    plot([-1.5,1.5],[-1.5,1.5],'-black');
    %axis([-1 1 -1 1 ]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    plot(c_stimulus(start_t + 30 - 1), c_stimulus(start_t + 29 - 1), 'ro','MarkerEdgeColor','red','MarkerSize',7,'LineWidth',2.5); % optimal stimulus
    plot(c_stimulus(start_t + 29 - 1:start_t + 30 - 1), [c_stimulus(start_t + 29 - 1),c_stimulus(start_t + 29 - 1)],'Color','red','LineWidth',2.5); % optimal stimulus
    
    
    plot(c_stimulus(start_t + 32 - 1), c_stimulus(start_t + 31 - 1),'ro','MarkerEdgeColor','black','MarkerSize',7,'LineWidth',2.5); % opposite stimulus
    plot(c_stimulus(start_t + 31 - 1:start_t + 32 - 1), [c_stimulus(start_t + 31 - 1), c_stimulus(start_t + 31 - 1)],'Color','black','LineWidth',2.5); % optimal stimulus

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    set(gca,'XTick',-0.8:0.2:0.8,'YTick',-0.8:0.2:0.8);
    axis([-0.75 0.75 -0.75 0.75 ]);
    xlabel(['Stimuli for t']);
    ylabel(['Stimuli for t-1']);



width = 4*3.8;     % Width in inches for subplot(4,3,time_point);
height = 4*1;    % Height in inches

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
 
set(gcf,'Position', myfiguresize);
set(gcf,'Resize','off');
% Save the file as PNG
print('Figure_Explanation','-dpng','-r300');
print('Figure_Explanation LOW_RES','-dpng','-r100');
end

%% Support functions
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