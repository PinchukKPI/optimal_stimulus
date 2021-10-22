function Fugire_optimal_stimulus
% Figure for the paper "Brightness change is optimal stimulus for 
% magnocellular-projecting retinal ganglion cells" by Artem Pinchuk

clc;
clear all;

alphabet = ['A','B','C','D','E','F','G','H','I','J','K' ];

figure('Name','Optimal stimulus','NumberTitle','off','Color', 'w', 'units','inches'); 
clf;

%% Delay responce for short stimulus
subplot(1,3,[1 2 3]); hold on;

stimulus_mask = 0.5;
stimulus_mask_size = 5;
step = 50;
shift = 25;
all_time = 6*step;

stimulus_line = zeros(6*step, 1);
stimulus_time = zeros(6*step, 1);
repeat_counter = 0;

time_stemp = 0;

i = 0;
while i < all_time
    i = i + 1;
    time_stemp = time_stemp + 1;
    stimulus_time(i) = time_stemp;
    if i == 5 * step - shift + repeat_counter
        for index = 1:1*stimulus_mask_size
            stimulus_line(i + index) = stimulus_mask +0.5;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
    elseif i == 5 * step - shift + repeat_counter + stimulus_mask_size - 1
        for index = 1:stimulus_mask_size
            stimulus_line(i + index) = - stimulus_mask +0.5;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
    elseif i == 4 * step - shift + repeat_counter
        for index = 1:stimulus_mask_size
            stimulus_line(i + index) = stimulus_mask +0.25;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
    elseif i == 4 * step - shift + repeat_counter + stimulus_mask_size - 1
        for index = 1:stimulus_mask_size
            stimulus_line(i + index) = - stimulus_mask +0.25;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
        
    elseif i == 3 * step - shift + repeat_counter
        for index = 1:stimulus_mask_size
            stimulus_line(i + index) = stimulus_mask;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
    elseif i == 3 * step - shift + repeat_counter + stimulus_mask_size - 1
        for index = 1:stimulus_mask_size
            stimulus_line(i + index) = - stimulus_mask;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
        
    elseif i == 2 * step - shift + repeat_counter
        for index = 1:stimulus_mask_size
            stimulus_line(i + index) = stimulus_mask-0.25;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
     elseif i == 2 * step - shift + repeat_counter + stimulus_mask_size - 1
        for index = 1:stimulus_mask_size
            stimulus_line(i + index) = - stimulus_mask - 0.25;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
       
    elseif i == step - shift + repeat_counter
        for index = 1:stimulus_mask_size
            stimulus_line(i + index) = stimulus_mask - 0.5;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
    elseif i == step - shift + repeat_counter + stimulus_mask_size - 1
        for index = 1:stimulus_mask_size
            stimulus_line(i + index) = - stimulus_mask - 0.5;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
    end 
end


up_shift = 3;
% plots for legend 
plot(stimulus_time(28:34), stimulus_line(28:34)+ up_shift, 'r','LineWidth',3); % example for optimal stimulus
plot(stimulus_time(77:82), stimulus_line(77:82)+ up_shift, 'black','LineWidth',3); %  
lh = legend('Optimal DS','Opposite DS', 'Location','NorthWest');
set(lh, 'Box', 'off');
set(lh, 'Color', 'none');



 plot(stimulus_time(10:5*step + 5), stimulus_line(10:5*step + 5) + up_shift, 'b','LineWidth',1.5);
 
 % optimal stimulus
plot(stimulus_time(28:34), stimulus_line(28:34)+ up_shift, 'r','LineWidth',3); % example for optimal stimulus
plot(stimulus_time(82:88), stimulus_line(82:88)+ up_shift, 'r','LineWidth',3); %  
plot(stimulus_time(136:142), stimulus_line(136:142)+ up_shift, 'r','LineWidth',3); %  
plot(stimulus_time(190:196), stimulus_line(190:196)+ up_shift, 'r','LineWidth',3); %  
plot(stimulus_time(244:250), stimulus_line(244:250)+ up_shift, 'r','LineWidth',3); %  

% opposite stimulus

plot(stimulus_time(77:82), stimulus_line(77:82)+ up_shift, 'black','LineWidth',3); %  
plot(stimulus_time(131:136), stimulus_line(131:136)+ up_shift, 'black','LineWidth',3); %  
plot(stimulus_time(185:190), stimulus_line(185:190)+ up_shift, 'black','LineWidth',3); %  
plot(stimulus_time(239:244), stimulus_line(239:244)+ up_shift, 'black','LineWidth',3); % 

stimulus_mask = 1.0; 
stimulus_mask_size = 5;
stimulus_line = zeros(6*step, 1);
stimulus_time = zeros(6*step, 1);
repeat_counter = 0;

time_stemp = 0;

i = 0;
while i < all_time
    i = i + 1;
    time_stemp = time_stemp + 1;
    stimulus_time(i) = time_stemp;
    if i == 5 * step - shift + repeat_counter
        for index = 1:1*stimulus_mask_size
            stimulus_line(i + index) = stimulus_mask;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + 1*stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;

    elseif i == 4 * step - shift + repeat_counter
        for index = 1:2*stimulus_mask_size
            stimulus_line(i + index) = stimulus_mask;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + 2*stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
        
    elseif i == 3 * step - shift + repeat_counter
        for index = 1:3*stimulus_mask_size
            stimulus_line(i + index) = stimulus_mask;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + 3*stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
        
    elseif i == 2 * step - shift + repeat_counter
        for index = 1:4*stimulus_mask_size
            stimulus_line(i + index) = stimulus_mask;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + 4*stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
        
    elseif i == step - shift + repeat_counter
        for index = 1:5*stimulus_mask_size
            stimulus_line(i + index) = stimulus_mask;
            stimulus_time(i + index) = time_stemp;
            time_stemp = time_stemp + 1;
        end
        i = i + 5*stimulus_mask_size;
        time_stemp = time_stemp - 2;
        repeat_counter = repeat_counter + 2;
    end 
end


plot(stimulus_time(10:5*step), stimulus_line(10:5*step), 'b','LineWidth',1.5);
% optimal stimulus
plot(stimulus_time(48:53), stimulus_line(48:53), 'r','LineWidth',3); % example for optimal stimulus
plot(stimulus_time(95:100), stimulus_line(95:100), 'r','LineWidth',3); %  
plot(stimulus_time(142:147), stimulus_line(142:147), 'r','LineWidth',3); %  
plot(stimulus_time(189:194), stimulus_line(189:194), 'r','LineWidth',3); %  
plot(stimulus_time(236:241), stimulus_line(236:241), 'r','LineWidth',3); %  

% opposite stimulus
plot(stimulus_time(23:28), stimulus_line(23:28), 'black','LineWidth',3); % example for opposite stimulus
plot(stimulus_time(75:80), stimulus_line(75:80), 'black','LineWidth',3); %  
plot(stimulus_time(127:132), stimulus_line(127:132), 'black','LineWidth',3); %  
plot(stimulus_time(179:184), stimulus_line(179:184), 'black','LineWidth',3); %  
plot(stimulus_time(231:236), stimulus_line(231:236), 'black','LineWidth',3); %  
set(gca,'YGrid','on');
set(gca,'YLim', [-0.5, 5]);
set(gca,'YTick',-0:3:3);
set(gca,'YTickLabel',{'0', '0' });

set(gca,'XLim', [0, 250]);
set(gca,'XTick',0:250:250);
set(gca,'XTickLabel',{'          Fast responce', 'Slow responce           ' });

ylabel('Stimulus');





hold off;

width = 4*3;     % Width in inches for subplot(4,3,time_point);
height = 4*1;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;  % 11     % Fontsize
lw = 1.5;      % LineWidth
msz = 9;       % MarkerSize
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
print('Figure_optimal_stimulus','-dpng','-r600');
print('Figure_optimal_stimulus LOW_RES','-dpng','-r100');
%print('Figure_STBM','-dpng','-r600')


end