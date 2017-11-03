% calculate drift via lin reg over 100 seconds
x = [];
y = [];
time = [];

for i = 0:500
    while true
        im = dir('live_images\*.tif');
        
        if length(im) >= 2
            break;
        end
    end
    
    data = processing(im, 0, 0, 0, 0);
    meanVector = data.mean;
    
    x = [x; meanVector(1)];
    y = [y; meanVector(2)];
    time = [time; str2num(im(2).name(1:end-4))];
    delete(fullfile('live_images\', im(2).name));
end

xFit = polyfit(time, x, 1);
yFit = polyfit(time, y, 1);

xDrift = xFit(1)
yDrift = yFit(1)


