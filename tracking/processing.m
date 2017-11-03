function results = processing(images, bool, iteration, xDrift, yDrift)
%% cycle through and process given images (written for 2 images)
pos=[];
Ithresh=500;
objectSize=25;
FOLDER = 'live_images\'
% FOLDER = 'images\'  % just for manual analysis

% repeats process for both images
for i=1:2
    I =imread(fullfile(FOLDER, images(i).name)); %im2uint8()
    %I=imadjust(I);
    I=bpass(I,1,objectSize,Ithresh);
    %%%imshow(I)
    pks=pkfnd(I,Ithresh,objectSize);
    out=cntrd(I,pks,objectSize);
    pos=[pos; out(:,1:2) , (i-1)*ones(length(out(:,1)),1)];
    
%             figure
%             colormap('gray');
%             imshow(I)
%             hold all
%             plot(out(:,1),out(:,2),'o')
end
%%
param.mem=1;
param.dim=2;
param.good=0;
param.quiet=1;
trackout = track(pos,4,param);

%%
n=1;
pos0=[];
pos1=[];
%loop through all tracks
for i=1:trackout(end,4)
    indx=find(trackout(:,4)==n);
    if length(indx)==2
        pos0=[pos0;trackout(indx(1),1:2)];
        pos1=[pos1;trackout(indx(2),1:2)];
    end
    n=n+1;
end
% images
% pos1
% pos0
% test = trackout(indx(1),1:2)

%% correction for microscope drift
pos1(:, 1) = pos1(:, 1) - (xDrift*iteration);
pos1(:, 2) = pos1(:, 2) - (yDrift*iteration);

dX = pos1-pos0;


%% save mat file for every nth image
if bool == 1
    save(strcat('data\', images(2).name(1:end-4), '.mat'), 'pos0', 'pos1');
end
%%

results.mean = [];
results.mean = [mean(dX(:, 1)), mean(dX(:, 2))];
results.sd = std(sqrt(dX(:, 1).^2 + dX(:, 2).^2));
results.disp = dX;
results.pos = pos0;

end







