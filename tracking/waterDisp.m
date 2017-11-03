clear all;
% batch_main;

folder = 'data\';
files = dir(fullfile(folder, '*.mat'));
pressure = [];
disp = [];
sd = [];

for i=1:length(files)
    D = load(strcat(folder, files(i).name));
    pressure = [pressure; str2num(files(i).name(1:4))];
    disp = [disp; sqrt(D.meanDisp(1)^2 + D.meanDisp(2)^2)];
    sd = [sd; D.sd]; 
end

dispMeter = disp .* (6.5e-6/20);
scatter(pressure, dispMeter, 'filled');
% errorbar(pressure, disp, sd, 'o');
hold on;
% plot(fitlm(pressure, disp));
title('Displacement as a function of pressure');
xlim([0 2700]);
xlabel('Pressure [mbar]');
ylabel('Displacement [m]');
