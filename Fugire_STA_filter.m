function Fugire_STA_filter
% Figure for the paper "Brightness change is optimal stimulus for 
% magnocellular-projecting retinal ganglion cells" by Artem Pinchuk
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

%   CELL 13 non compressed AREA for center   above axis (113301819.2433)  and under axis (-112265314.6404) above/under = -1.0092
cellInd = 13;  % 1-8 Off-center   9-15 On-center

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
STA_time_period = 2000; % 2000 = 200 ms because sampleRate = 1e4; 10000 Hz  

compression_coef = 100;
c_stimulus = compress_stimuly(stimulus, compression_coef);
c_spike_response = compress_spikes(spike_response, compression_coef);
STA = STA_calc(c_stimulus, c_spike_response, floor(STA_time_period/compression_coef));
STA_sum = 0;
for i = 1:floor(STA_time_period/compression_coef)
    STA_sum = STA_sum + STA(i);
end
% figure(3);
% hold on
% plot(STA,'r');
% grid on;
% title(['SUM STA(',num2str(floor(STA_time_period/compression_coef)), ')=',num2str(STA_sum) ]);
% hold off
%% STA distribution
%  
compression_coef = 100;
STA_time_period = 2000;
c_stimulus = compress_stimuly(stimulus, compression_coef);
c_spike_response = compress_spikes(spike_response, compression_coef);
ho_size = size(c_spike_response);
new_size = ho_size(1); %    1195200 
STA_c_time_period = floor(STA_time_period / compression_coef); 
STA_c = STA_calc(c_stimulus, c_spike_response, STA_c_time_period);
%figure('Name','Spike-triggered average','NumberTitle','off','Color', 'w');
%hold on
% reverse_i = zeros(STA_c_time_period, 1);
% for i = -1:-1:-STA_c_time_period
%     reverse_i(abs(i)) = i;
% end


%plot(reverse_i, STA_c, '--rs','LineWidth',2);


%grid(gca,'minor');
%hold off;

%  CELL 13 AREA for center   above axis (113301819.2433)  and under axis (-112265314.6404) -1.0092
% after compression 
c_area_above = 0;
c_area_under = 0;
summ = 0;
for i = 1:STA_c_time_period
    summ = summ + STA_c(i);
    if STA_c(i) > 0
        c_area_above = c_area_above + STA_c(i);
    else
        c_area_under = c_area_under + STA_c(i);
    end
end
disp(['cell #',num2str(cellInd), ' compressed AREA for center   above axis (',num2str(c_area_above), ')  and under axis (',num2str(c_area_under), ') ', num2str(c_area_above/c_area_under) ]);

for i = 1:STA_c_time_period
    STA_c(i) = STA_c(i) - summ / STA_c_time_period;  %normalization trick to look more like derivative kernel 
    % because after compression areas more different 
    % non compressed areas for CELL 13 AREA for center   above axis (113301819.2433)  and under axis (-112265314.6404) above/under =-1.0092
end

c_area_above = 0;
c_area_under = 0;
for i = 1:STA_c_time_period
    if STA_c(i) > 0
        c_area_above = c_area_above + STA_c(i);
    else
        c_area_under = c_area_under + STA_c(i);
    end
end
disp(['cell #',num2str(cellInd), ' compressed and norm AREA for center   above axis (',num2str(c_area_above), ')  and under axis (',num2str(c_area_under), ') ', num2str(c_area_above/c_area_under) ]);

disp(' ');
%% stimulus and STA as filtet

    % CenterSurroundWhiteNoise is a cell array, each entry of which is a structure that corresponds
    % to a cell in the dataset. Structure fields:
    % .stimulus = concatenated stimulus traces for .center and .surround
    % stimuli
    % .response = concatenated excitatory conductance response traces (in nS)
    % for .center, .surround and .centerSurround stimuli
    % Note that data are concatenated and grouped for convenience, but were
    % acquired in interleaved trials 

    on_correlation_center_above = zeros(7,1); 
    on_correlation_surround_above = zeros(7,1); 
    on_correlation_center_under = zeros(7,1); 
    on_correlation_surround_under = zeros(7,1); 
    off_correlation_center_above = zeros(8,1); 
    off_correlation_surround_above = zeros(8,1); 
    off_correlation_center_under = zeros(8,1); 
    off_correlation_surround_under = zeros(8,1); 
    
    alphabet = ['A','B','C','D','E','F','G','H','I','J','K' ];
    empty_space   = '                                                                               ';
    empty_space_2 = empty_space;
    empty_space_3 = empty_space;
    gray_color = [.45 0.45 0.45];
    %figure('Name','STA as filtet','NumberTitle','off','Color', 'w', 'units','normalized','outerposition', [0 0 1 1]); 
    figure('Name','STA as filtet','NumberTitle','off','Color', 'w', 'units','inches'); 
    clf;
    set(gcf,'Resize','off');
    
    for cc = 1:15
        cellInd = cc;
        sampleRate = 1e4; %Hz 10000
        filterLen = 500; %msec, length of linear filter to compute

        %frequency after which to cut off filter spectrum. The stimulus updated at 30
        %Hz, so the cutoff should be below that value to avoid ringing
        
        %using code of function getLinearFilter() from paper
        % "Receptive field center-surround interactions mediate context-dependent 
        % spatial contrast encoding in the retina"
        % Maxwell H Turner, Gregory W Schwartz, Fred Rieke
        % DOI: https://doi.org/10.7554/eLife.38841 
        % https://elifesciences.org/articles/38841
        filterPts = (filterLen/1000)*sampleRate; %msec -> datapoints

        centerFilter = getLinearFilter(CenterSurroundWhiteNoise{cellInd}.stimulus.center,...
            CenterSurroundWhiteNoise{cellInd}.response.center);
        centerFilter = centerFilter(1:filterPts); %trim to just filterLen
        %centerFilter = fliplr(centerFilter);
        surroundFilter = getLinearFilter(CenterSurroundWhiteNoise{cellInd}.stimulus.surround,...
            CenterSurroundWhiteNoise{cellInd}.response.surround);
        surroundFilter = surroundFilter(1:filterPts);  %trim to just filterLen
        %surroundFilter = fliplr(surroundFilter);
        filterTimeVector = -(1:filterPts) ./ sampleRate * 1000; %sec
        %filterTimeVector = -(filterPts:-1:1) ./ sampleRate; %sec

        
        
        if cc < 9
            subplot(3,3,1, 'YAxisLocation', 'right'); hold on;
            plot(filterTimeVector, centerFilter,'b', 'LineWidth',1.5);
            plot(filterTimeVector, surroundFilter,'r', 'LineWidth',1.5);
        else
            subplot(3,3,4, 'YAxisLocation', 'right'); hold on;
            plot(filterTimeVector, centerFilter,'Color','m', 'LineWidth',1.5);  % [.6 0 .6]
            plot(filterTimeVector, surroundFilter,'Color',gray_color,'LineWidth',1.5);  %[.65 0.65 0.65]
        end
        
        c_area_above = 0;
        c_area_under = 0;
        s_area_above = 0;
        s_area_under = 0;
        for i = 1:filterPts
            if centerFilter(i) > 0
                c_area_above = c_area_above + centerFilter(i);
            else
                c_area_under = c_area_under + centerFilter(i);
            end
            if surroundFilter(i) > 0
                s_area_above = s_area_above + surroundFilter(i);
            else
                s_area_under = s_area_under + surroundFilter(i);
            end
        end
        if cc < 9
            off_correlation_center_above(cc) = c_area_above; 
            off_correlation_surround_above(cc) = s_area_above; 
            off_correlation_center_under(cc) = c_area_under; 
            off_correlation_surround_under(cc) = s_area_under; 
        else
            on_correlation_center_above(cc-8) = c_area_above; 
            on_correlation_surround_above(cc-8) = s_area_above; 
            on_correlation_center_under(cc-8) = c_area_under; 
            on_correlation_surround_under(cc-8) = s_area_under; 
        end
        disp(['  CELL ',num2str(cc) ,' AREA center   above axis (',num2str(c_area_above), ') under axis (',num2str(c_area_under), ') above/under =', num2str(c_area_above/c_area_under) ]);
        disp(['  CELL ',num2str(cc) ,' AREA surround above axis (',num2str(s_area_above), ') under axis (',num2str(s_area_under), ') above/under =', num2str(s_area_above/s_area_under) ]);
        disp(' ');

    end
    
    
    zero_line = zeros(filterPts, 1); 
    
    hold off
	subplot(3,3,1); hold on;
    plot(filterTimeVector, zero_line,'Color','black','LineWidth',1);
	title({['\bf\fontsize{11}', alphabet(1), empty_space];['\rm\fontsize{10} OFF ']});
    legend('Off Center','On Surround', 'Location','SouthWest');
    xlabel('Time [ms]');
    ylabel('Temporal receptive field');
    %set(gca,'YTick',0);
    set(gca,'YTickLabel',{' '});
    lh = legend();
    set(lh, 'Box', 'off');
    set(lh, 'Color', 'none');
    
    subplot(3,3,4); hold on;
    plot(filterTimeVector, zero_line,'Color','black','LineWidth',1);
 	title({['\bf\fontsize{11}', alphabet(2), empty_space];['\rm\fontsize{10} ON ']});
    legend('On Center','Off Surround', 'Location','NorthWest');
    xlabel('Time [ms]');
    ylabel('Temporal receptive field');
    %set(gca,'YTick',0);
    % delete Exponent Label Using Ruler Objects 
    set(gca,'YTickLabel',{' '});
    lh = legend();
    set(lh, 'Box', 'off');
    set(lh, 'Color', 'none');
    
    subplot(3,3,[2 3 5 6]); hold on;
    hold on

    axis([0 400000000 0 400000000 ]);
    plot(abs(off_correlation_center_under),abs(off_correlation_center_above),'+','MarkerEdgeColor','b','MarkerSize',10, 'LineWidth',2);
    plot(abs(off_correlation_surround_under),abs(off_correlation_surround_above),'o','MarkerEdgeColor','r','MarkerSize',10, 'LineWidth',2);
    plot(abs(on_correlation_center_under),abs(on_correlation_center_above),'x','MarkerEdgeColor','m','MarkerSize',10, 'LineWidth',2);
    plot(abs(on_correlation_surround_under),abs(on_correlation_surround_above),'s','MarkerEdgeColor',gray_color,'MarkerSize',10, 'LineWidth',2);
 	title({['\bf\fontsize{11}', alphabet(3), empty_space_2,empty_space_2];['\rm\fontsize{10} Temporal receptive field area comparision ']});
    lh = legend('Off Center','Off Surround','On Center','On Surround', 'Location','NorthWest');
    set(lh, 'Box', 'off');
    set(lh, 'Color', 'none');
    plot([0,400000000],[0,400000000],'-red');
    xlabel('Area under axis');
    ylabel('Area above axis');
    set(gca,'YTickLabel',{'0', ' ', '2', ' ', '4', ' ', '6', ' ', '8' });
    set(gca,'XTickLabel',{' ', ' ', '2', ' ', '4', ' ', '6', ' ', '8' });
    hold off
% stimulus_mask = zeros(STA_c_time_period, 1);


step = 60;
stimulus_time_line = zeros(6*step, 1);
for i = 1:6*step
    if i > 4 * step
        stimulus_time_line(i) = 3.0;
    elseif i > 3 * step
        stimulus_time_line(i) = 1.5;
    elseif i > 2 * step
        stimulus_time_line(i) = 0.0;
    elseif i > step
        stimulus_time_line(i) = -1.5;
    elseif i > 0
        stimulus_time_line(i) = -3.0;
    end 
end

new_size = 6*step;

 for i = 1:new_size
    kernel_sum = 0;
    for filter_index = 1:STA_c_time_period
        if (i + filter_index - 24) > 0 && (i + filter_index - 1) < new_size
            kernel_sum = kernel_sum + STA_c(STA_c_time_period - filter_index + 1)* stimulus_time_line(i + filter_index - 24);  % - 12
        end
    end
    filtered(i) = kernel_sum;
 end

time_line = (1:new_size)*0.01;
 
subplot(3,3,[7 8 9]); hold on;
title({['\bf\fontsize{11}', alphabet(4), empty_space_3, empty_space_3,empty_space_3];['\rm\fontsize{10} ON cell center STA as filter ']});
plot(time_line(24:5*step), stimulus_time_line(24:5*step), 'r','LineWidth',1.5);

ax1 = gca;
set(ax1,'XLim', [0.25, 2.75]);

set(ax1,'YGrid','on');
set(ax1,'YLim', [-4, 4]);
set(ax1,'YTick',-4:2:4);
plot(time_line(24:5*step), filtered(24:5*step)*2,'b','LineWidth',1.5);
lh = legend('Stimulus','Responce', 'Location','NorthWest');
set(lh, 'Box', 'off');
set(lh, 'Color', 'none');
ylabel('Stimulus');
xlabel('Time [s]');

ax2 = axes('Position',get(ax1,'Position'),...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YLim', [-0.5, 0.5]);
set(ax2,'YTick',-0.5:0.25:0.5);
set(ax2,'XTick',[]);
set(get(ax2,'YLabel'),'String','Spike rate')



hold off;

% Defaults for this blog post
width = 4*3;     % Width in inches for subplot(4,3,time_point);
height = 4*3;    % Height in inches

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
print('Figure_STA_filter','-dpng','-r300');
print('Figure_STA_filter LOW_RES','-dpng','-r100');
%print('Figure_STBM','-dpng','-r600');

end

%%  Support functions

function LinearFilter = getLinearFilter(stimulus, response)
% Quickly estimate a linear filter using a white noise stimulus and response
FilterFft = mean((fft(response,[],2).*conj(fft(stimulus,[],2))),1) ;
LinearFilter = real(ifft(FilterFft)); %fliplr

end

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