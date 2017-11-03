function results = processing0(images)
%% cycle through and process given images (written for 2 images)

pos=[];
Ithresh=500;
objectSize=25;
FOLDER = 'live_images\'

% processes 0 mbar image

I =imread(fullfile(FOLDER, images(1).name)); %im2uint8()
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

%%
param.mem=1;
param.dim=2;
param.good=0;
param.quiet=1;
trackout = track(pos,4,param);
n=1;
%%

pos0=[];
%loop through all tracks
for j=1:trackout(end,4)
    indx=find(trackout(:,4)==n);
    if length(indx)==2
        pos0=[pos0;trackout(indx(1),1:2)];
    end
    n=n+1;
end
% images
% pos1
% pos0
% test = trackout(indx(1),1:2)
results = pos0;








