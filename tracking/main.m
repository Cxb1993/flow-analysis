names={'280 mbar off.2_MMStack.ome.tif','280 mbar on.2_MMStack.ome.tif'};
pos=[];
Ithresh=500;
objectSize=25;
i=1;
for i=1:length(names)
    I =(imread(fullfile('images\', names{i})));%im2uint8()
    %I=imadjust(I);
    
    I=bpass(I,1,objectSize,Ithresh);
    imshow(I)
    pks=pkfnd(I,Ithresh,objectSize);
    out=cntrd(I,pks,objectSize);
    pos=[pos; out(:,1:2) , (i-1)*ones(length(out(:,1)),1)];
    
    figure
    colormap('gray');
    imshow(I)
    hold all
    plot(out(:,1),out(:,2),'o')
end
%%
param.mem=1;
param.dim=2;
param.good=0;
param.quiet=1;
trackout = track(pos,4,param)
n=1;
%%
pos0=[];
pos1=[];
%loop through all tracks
for i=1:trackout(end,4)
    indx=find(trackout(:,4)==n);
    if length(indx)==2
        pos0=[pos0;trackout(indx(1),1:2)];
        pos1=[pos1;trackout(indx(2),1:2)];
    end
    n=n+1
end
pos1
pos0
dX=pos1-pos0;


imshow(I)
hold all
quiver(pos0(:,1),pos0(:,2),dX(:,1),dX(:,2))
%%
% save(strcat(names{1}(1:end-4),'.mat'),'pos0','pos1')
save(strcat('data\', names{1}(1:end-4), '.mat'), 'pos0', 'pos1');
scatter(dX(:,1), dX(:,2))




