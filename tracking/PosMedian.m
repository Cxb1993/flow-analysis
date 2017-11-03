function results = PosMedian(numberOfPoints)
% THIS PROGRAM IS CRAP
% function takes in the number of points from which to calculate the median
% position of all the particles

n = numberOfPoints;

while true
    im = dir('live_images\*.tif');
    
    if length(im) >= n
        break;
    end
end

Ithresh=500;
objectSize=25;
FOLDER = 'live_images\'

%% processes each image
for i = 1:n
    I =imread(fullfile(FOLDER, im(i).name)); %im2uint8()
    I=bpass(I,1,objectSize,Ithresh);
    pks=pkfnd(I,Ithresh,objectSize);
    out=cntrd(I,pks,objectSize);
    
    % out(:, 1:2) spits out x and y coordinates for each particle found
    % 0*ones(...) assigns a value of 0 to indicate these points are the
    % initial points (compatability with processing(images, bool))   
    calcPos = [out(:,1:2), 0*ones(length(out(:,1)),1)];
    
    coord(i).x = calcPos(:, 1);
    coord(i).y = calcPos(:, 2);
end

length(coord(1).x)
length(coord(7).x)

allX = cat(3, coord(1).x, coord(2).x, coord(3).x, coord(4).x, coord(5).x, coord(6).x, coord(7).x, coord(8).x, coord(9).x, coord(10).x);
allY = cat(3, coord(1).y, coord(2).y, coord(3).y, coord(4).y, coord(5).y, coord(6).y, coord(7).y, coord(8).y, coord(9).y, coord(10).y);

finalX = median(allX, 3);
finalY = median(allY, 3);


%% creates a temporary array to store similar coords to determine median

end





