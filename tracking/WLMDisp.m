clear all;
batch_main;

folder = 'data\';
files = dir(fullfile(folder, '*.mat'));
pressure = [];
disp = [];
direction = [];
sd = [];

for i=1:length(files)
    D = load(strcat(folder, files(i).name));
    pressure = [pressure; str2num(files(i).name(1:4))];
    disp = [disp; sqrt(D.meanDisp(1)^2 + D.meanDisp(2)^2)];
    direction = [direction; D.meanDisp];
    sd = [sd; D.sd]; 
end

% scatter(pressure, disp);
errorbar(pressure, disp, sd, 'o');
hold on;
% plot(fitlm(pressure, disp));
title('Displacement as a function of pressure with standard deviation bars');
xlim([250 300]);
xlabel('Pressure [mbar]');
ylabel('Displacement [px]');

% graph of directions of displacement
quiver(pressure, disp, direction(:, 1), direction(:, 2), 'r');



%% graphs stress
cd ..;

for i=1:length(files)  
    zeroPressure = 'tracking\images\0000 mbar 1.0_MMStack_Pos0.ome.tif';
    master(files(i).name, zeroPressure);
    % need to add automation of PIV processing
    hold on;
end