NUM = 4000;         % number of images to take

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;

xDrift = -5.5460e-04;
yDrift = -1.9631e-04;

for j = 1:NUM
    % for i = 1:NUM
    % naming convention:
    % 0 mbar image must start with '0's
    % all other images must have nonzero names
    while true
        im = dir('live_images\*.tif');
        
        if length(im) >= 2
            break;
        end
    end
    
    %im1 is second image
    im1 = im(2).name(1:end-4);
    i = str2num(im1);
    
    %% processes two images and extracts mean displacement
    n = 10;
    
    if mod(i, n) == 1
        % plot stress distribution for every nth graph
        % pass 1 for true and initiate stress graph
        % NOTE: you have to manually determine the x and y drift via linear
        % regression (should do over 100s)  
        % processing(images, bool, xDrift, yDrift)
        data = processing(im, 1, j, xDrift, yDrift);
        
        cd ..;
        iteration = num2str(i);      % converts i to number
        figure;
%         master(strcat('tracking\data\', iteration, '.mat'), strcat('tracking\live_images\', iteration, '.tif'), iteration);
        master(strcat('tracking\data\', iteration, '.mat'), strcat('tracking\live_images\0.tif'), iteration);
        xlabel('[px]');
        ylabel('[px]');
        drawnow;
        cd tracking;
        movefile(strcat('live_images\', im(2).name), 'saved_images');
    else
        data = processing(im, 0, j, xDrift, yDrift);
        delete(fullfile('live_images\', im(2).name));
    end
    
    meanVector = data.mean;
    meanDisp = sqrt(meanVector(1)^2 + meanVector(2)^2);
%     sd = data.sd;

     %% scatter plot
%     set(0, 'CurrentFigure', f1);
%     hold on;
%     scatter(i, meanDisp, 'o', 'b');
%     xlabel('Time [s]');
%     ylabel('Displacement [px]');
%     title('Mean displacement as a function of time (260 mbar)');
% %     hold on;
% %     quiver(i, meanDisp, meanVector(1), meanVector(2));
%     drawnow;
%     %quiver(t, meanDisp, meanVector(1), meanVector(2), 'b');
     %% scatter plot of x and y
    set(0, 'CurrentFigure', f1);
    hold on;
    scatter(i, meanVector(1), 'o', 'b');
    xlabel('Time [s]');
    ylabel('Displacement [px]');
    title('x Mean displacement as a function of time (260 mbar)');
    drawnow;
    
    set(0, 'CurrentFigure', f2);
    hold on;
    scatter(i, meanVector(2), 'o', 'b');
    xlabel('Time [s]');
    ylabel('Displacement [px]');
    title('y Mean displacement as a function of time (260 mbar)');
    drawnow;
    
    set(0, 'CurrentFigure', f4);
    hold on;
    scatter(i, meanDisp, 'o', 'b');
    xlabel('Time [s]');
    ylabel('Displacement [px]');
    title('Mean displacement as a function of time (260 mbar)');
    drawnow;

    
    %% angle over time
    set(0, 'CurrentFigure', f3);
    hold on;
    angle = atan2(meanVector(2), meanVector(1));
    scatter(i, angle, 'o', 'b');
    xlabel('Time[s]');
    ylabel('Angle [rad]');
    title('Angle with respect to time');

end
