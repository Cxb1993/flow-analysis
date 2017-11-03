IMAGE1 = '0 mbar 20x_MMStack_Pos0.ome.tif';
IMAGE2 = '260 mbar 20x_MMStack_Pos0.ome.tif';


% processing for IMAGE1
a1 = double(imread(IMAGE1));
b1 = bpass(a1, 1, 20);

% locates particles for IMAGE1
temp1 = pkfnd(b1, 25, 8);
temp1(:, 3) = 0;
num = numel(temp1);

% processing for IMAGE2
a2 = double(imread(IMAGE2));
b2 = bpass(a2, 1, 20);

% locates particles for IMAGE2
temp2 = pkfnd(b2, 25, 8);
temp2(:, 3) = 1;
position = [temp1; temp2];

%% centroid calculation
centroids1 = cntrd(b1, temp1, 7);
centroids2 = cntrd(b2, temp2, 7);

% clunky but works
centroids1_copy = centroids1(:, 1:2);
centroids1_copy(:, 3) = 0;
centroids2_copy = centroids2(:, 1:2);
centroids2_copy(:, 3) = 1;
position_detailed = [centroids1_copy; centroids2_copy];

%% average distance
% lastPoint = length(temp2);
% totalDistance = 0;
% for i = 1:(lastPoint-1)
%     for j = i+1:lastPoint
%         difference = temp2(i, 1:2)-temp2(j, 1:2);
%         distance = sqrt(difference(1)^2 + difference(2)^2);
%         totalDistance = totalDistance + distance;
%     end
% end
% 
% avgDistance = totalDistance/(0.5*(lastPoint-1)^2);

%% tracking
MAX_DISPLACEMENT = 10;
delete = [];
tracks = track(position, MAX_DISPLACEMENT);
tracks_advanced = track(position_detailed, MAX_DISPLACEMENT);

for i = 1:(size(tracks, 1)-1)
    if tracks(i, 3) == tracks(i+1, 3) && tracks(i, 3) == 0
        delete = [delete, i];
    elseif tracks(i, 3) == tracks(i+1, 3) && tracks(i, 3) == 1
        delete = [delete, i+1];
    end
end

for i = 1:(size(tracks_advanced, 1)-1)
    if tracks_advanced(i, 3) == tracks_advanced(i+1, 3) && tracks_advanced(i, 3) == 0
        delete = [delete, i];
    elseif tracks_advanced(i, 3) == tracks_advanced(i+1, 3) && tracks_advanced(i, 3) == 1
        delete = [delete, i+1];
    end
end

tracks(delete, :) = [];
tracks_advanced(delete, :) = [];
 
%% tracking visualization
coords = [];
coords_shifted = [];

for i = 1:size(tracks_advanced, 1)
    if tracks_advanced(i, 3) == 0
        coords = cat(1, coords, tracks_advanced(i, :));
    else 
        coords_shifted = cat(1, coords_shifted, tracks_advanced(i, :));
    end
end

d(1).r = coords;
d(2).r = coords_shifted;
%% creating a grid

% Subtract off displacements from reference time
displacement = d(2).r-d(1).r;

% Select number of points for interpolated grid
ovr = 1; % Spatial oversampling (ovr=1 gives grid spacing= avg interparticle distance). ovr should be <=1.
beadNumber = length(d(1).r); % Total number of beads
nx = round(ovr*sqrt(beadNumber)); % Number of points on each side of the interpolation grid
if mod(nx, 2) == 0
    nx = nx+1; % Make sure odd number points in grid
end

% We pad the displacement data on each side of the grid. This reduces
% artifacts in the stress calculation. fracpad is the fraction of extra padding on
% each side of the original data. Thus fracpad=0.5 doubles the width and
% height of the original field of view
fracpad=0.5;
npad = round(fracpad*nx);

% Calculate the boundaries of the data set
xmn = min(d(1).r(:,1));
xmx = max(d(1).r(:,1));
ymn = min(d(1).r(:,2));
ymx = max(d(1).r(:,2));

dx = max((xmx-xmn)/nx, (ymx-ymn)/nx); % Distance between the grid points
c = 0.5*[xmn+xmx, ymn+ymx]; % Centre of data set

% Construct the grid
xi = linspace(-(nx-1)/2-npad,(nx-1)/2+npad,nx+2*npad)*dx+c(1);
yi = linspace(-(nx-1)/2-npad,(nx-1)/2+npad,nx+2*npad)*dx+c(2);
[X,Y] = meshgrid(xi,yi); % Matrix of gridpoints

%% Interpolate the particle track data onto the grid
% If adapting to 3d TFM, the out-of-plane displacements should be
% interpolated onto the X,Y grid to give d(i).dz_interp.
%pix = 1e-6; % Size of one pixel in meters
NUM_PARTICLES = 8;

% surface_interpolate(x,y,z,xi,yi,num)
% the bigger num is, the more smoothing that is applied; excess smoothing
% can be avoided by adjusting size of particles that pkfnd searches for
d(2).dx_interp = surface_interpolate(d(2).r(:,1), d(2).r(:,2), displacement(:,1), X, Y, NUM_PARTICLES);
d(2).dy_interp = surface_interpolate(d(2).r(:,1), d(2).r(:,2), displacement(:,2), X, Y, NUM_PARTICLES);
%d(i).dz_interp=surface_interpolate(d(i).r(:,1),d(i).r(:,2),d(i).dr(:,3),X,Y,NUM_PARTICLES);
%find indices of all NaNs
ind2 = find(isnan(d(2).dx_interp));

% Calculate field of view array at each timepoint (points within the
% original data set). fov has zeros outside of field of view, and ones
% inside.
d(2).fov = ones(size(X));
d(2).fov(ind2) = 0;

imshow(IMAGE1);
hold all;
if i>1
    %quiver(X/pix,Y/pix,d(2).dx_interp/pix,d(2).dy_interp/pix);
    SCALE = 5;
    quiver(X, Y, d(2).dx_interp, d(2).dy_interp, SCALE);
    axis image
    xlabel('x [pixels]');
    ylabel('y [pixels]');
    title('Map of interpolated displacements')
    %axis([0 size_cell_image(2) 0 size_cell_image(1)]);
    pause(0.1)
end

%contourf(X/pix,Y/pix,d(i).dz_interp/pix,10)
%axis image

%% quiver plot
xShift = displacement(:, 1);
yShift = displacement(:, 2);
quiver(d(1).r(:, 1), d(1).r(:, 2), xShift, yShift, SCALE);

%% visualization
%  colormap('gray');
%  image(b1);
%  hold on;
%  pk = pkfnd(b1,5,5);
%  plot(pk(:,1), pk(:,2),'rd');