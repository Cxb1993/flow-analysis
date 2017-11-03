function master(dataName, imageName, iteration)
% Take particle tracks from images, and interpolate the particle tracks to
% get displacements on a regular grid, and finally calculate stresses.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                                                         %
% %      Example code for performing TFM. Takes particle tracks from images %
% %      and calculates stresses. Compiled by Rob Style, 21st Jan 2014      %
% %                                                                         %
% %      This code is for 2d traction force microscopy - generally suitable %
% %      for cellular TFM. This assumes that there are no out-of-plane      %
% %      tractions exerted on the surface of the substrate. In this case,   %
% %      only in-plane displacements at the surface of the substrate are    %
% %      required. The code can be easily modified for 3d TFM. The main     %
% %      important difference is in calculating the Q matrix, as noted in   %
% %      the code.                                                          %
% %                                                                         %
% %      This is based on various codes written by members of the Soft      %
% %      Matter Lab at Yale including Eric Dufresne, Ye Xu                  %
% %      Callen Hyland, Aaron Mertz and Ross Boltyanskiy                    %
% %                                                                         %
% %      Experimental data provided by Callen Hyland. Shows an Aplysia bag  %
% %      cell neuron growth cone on a soft, silicone gel substrate. The     %
% %      scale is 0.116 microns per pixel.                                  %
% %                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% example_tracks contains a structure array d that contains the tracked
% particle data from (i) a reference image, and (ii) an image where tractions are
% being applied to the surface of the TFM substrates.
% In this example, image1.tif is the reference image of fluorescent beads
% on a bare substrate. image2.tif shows the fluorescent beads with a aplysia bag cell
% growth cone sitting on the substrate.
% d(1).r gives the positions of the fluorescent beads in image1.tif. d(2).r
% gives the positions of the fluorescent beads in image2.tif with drift
% subtracted. d(1).dr is the displacements of the fluorescent beads in
% image1.tif relative to those in image1.tif (i.e. zero for all beads).
% d(2).dr gives the displacements of the fluorescent beads in image2.tif
% relative to those in image1.tif, with drift subtracted.
% In general, this TFM code requires a structure array input d. Each structure
% d(i) in the structure array should contain the tracked particle data for
% a specific time point. d(i).r should be the positions of all the tracked
% particles at time i. d(i).dr should be the displacement of these
% particles from their position at the first time point. Note that only
% data for particles that are tracked through every time point should be
% included.
%
% *Note* This code assumes that the fluorescent beads are at the surface
% of the substrate - not below the surface. We do not recommend that beads
% at a fixed distance below the surface are used to calculate surface
% traction stresses (see manuscript). However, if this is done, the Q
% matrix calculation below needs to be suitably modified by changing the
% inputs to calcQ. If strain energy density is required, it must be
% calculated with values of displacements at the surface - not below the
% surface. These surface displacements can be estimated from sub-surface
% displacements, as described in Mertz et al., PRL (2012) 108, 198101.
%
% The code has various outputs:
% X, Y are the x, y coordinates of the grid that stresses and displacements
% are returned at.
%
% d(i).dx_interp, d(i).dy_interp are the displacements d(i).dr interpolated
% onto the grid given by X,Y. Given at each time i.
%
% d(i).ux, d(i).uy are the interpolated displacements above, with a low-pass
% exponential filter applied. The filter is necessary to reduce noise, and
% is controlled by min_feature_size, described below. Given at each time i.
%
%
% d(i).fov is the field of view. Has the same size as X,Y. This is 1 when
% inside the area where we perform particle tracking, and is 0 outside
% (e.g. in the extra padded region that we add in the code). Given at each
% time i.
%
% d(i).stress_x, d(i).stress_y are the surface traction stresses calculated
% from the displacements d.dr. These are returned at the gridpoints given
% by X,Y. Given at each time i.
%
% d(i).sed is the strain energy density - the work per unit area required
% to deform the substrate by the cell/object of interest. See Mertz et al.,
% PRL 108, 198101 (2012) for details.

%load ('cell_outline.mat')

% cell_outline contains a series of points on the periphery of the cell.
% These are clicked by hand. plot(cell_x,cell_y) will show the cell outline
% in pixel units.

tref=1;  % Time for 'zero stress' reference. e.g. if a cell/object of
% interest is removed from the surface at time point 1, then tref=1. If the
% cell/object of interest is present and removed in the last frame,
% tref=length(d).

% This is the spatial resolution of the stress measurement in units of the grid spacing.
% The smaller min_feature_size the better your spatial resoltuion, but
% the worse signal to noise in the stress.
% System parameters
pix = 6.5e-6/20; %13e-6./20;  Size of one pixel in meters
EM = 3e+3; % Young's modulus in Pascal for Dow Corning CY 52-276 A/B
h = 45e-6; %37e-6;  % Film thickness in microns
nu =.499; % Poisson's ratio. If Poisson's ratio=1/2, set nu=0.499 to avoid some division by (1-2*nu) issues in the code.
DATA = dataName;
IMAGE = imageName;
%PIV = 'PIV_0705.mat';

%P=load(strcat('tracking\data\', DATA));
P = load(DATA);
%Inputs for displacements
%Z0=h*ones(1000,1);
d(1).r=P.pos0*pix;
d(2).r=P.pos1*pix;
d(1).dr=P.pos0*pix*0;
d(2).dr=d(2).r-d(1).r;
%% Create the grid to interpolate the particle tracking data onto
% This section should remain unchanged if adapting to 3d TFM.

% Subtract off displacements from reference time
for i=1:length(d)
    d(i).dr=d(i).dr-d(tref).dr;
end

% Select number of points for interpolated grid
ovr = 1; % Spatial oversampling (ovr=1 gives grid spacing= avg interparticle distance). ovr should be <=1.
nb_beads=length(d(1).r); % Total number of beads
nx=round(ovr*sqrt(nb_beads)); % Number of points on each side of the interpolation grid
if mod(nx,2)==0
    nx = nx+1; % Make sure odd number points in grid
end

% We pad the displacement data on each side of the grid. This reduces
% artifacts in the stress calculation. fracpad is the fraction of extra padding on
% each side of the original data. Thus fracpad=0.5 doubles the width and
% height of the original field of view
fracpad=0.5;
npad = round(fracpad*nx);
% Calculate the boundaries of the data set
xmn = min(d(tref).r(:,1));
xmx = max(d(tref).r(:,1));
ymn = min(d(tref).r(:,2));
ymx = max(d(tref).r(:,2));

dx = max( (xmx-xmn)/nx, (ymx-ymn)/nx); % Distance between the grid points
c=.5*[xmn+xmx,ymn+ymx]; % Centre of data set

% Construct the grid
xi = linspace(-(nx-1)/2-npad,(nx-1)/2+npad,nx+2*npad)*dx+c(1);
yi = linspace(-(nx-1)/2-npad,(nx-1)/2+npad,nx+2*npad)*dx+c(2);
[X,Y]=meshgrid(xi,yi); % Matrix of gridpoints

%% Interpolate the particle track data onto the grid
% If adapting to 3d TFM, the out-of-plane displacements should be
% interpolated onto the X,Y grid to give d(i).dz_interp.

%figure
im = double(imread(IMAGE));
%colormap gray
%hold on
size_cell_image=size(im); % Size of control image in pixels

%imagesc(im);
% size(d(1).dr(:,1));
% size(d(1).dr(:,2));
%size(d(1).dr(:,3))
for i = 1:length(d)
    
    d(i).dx_interp=surface_interpolate(d(i).r(:,1),d(i).r(:,2),d(i).dr(:,1),X,Y,8);
    d(i).dy_interp=surface_interpolate(d(i).r(:,1),d(i).r(:,2),d(i).dr(:,2),X,Y,8);
    %d(i).dz_interp=surface_interpolate(d(i).r(:,1),d(i).r(:,2),d(i).dr(:,3),X,Y,8);
    %find indices of all NaNs
    ind2=find(isnan(d(i).dx_interp));
    
    % Calculate field of view array at each timepoint (points within the
    % original data set). fov has zeros outside of field of view, and ones
    % inside.
    d(i).fov=ones(size(X));
    d(i).fov(ind2)=0;
    
    %     if i>1
    %         quiver(X/pix,Y/pix,d(i).dx_interp/pix,d(i).dy_interp/pix);
    %         axis image
    %         xlabel('x [pixels]');
    %         ylabel('y [pixels]');
    % %         title('Map of interpolated displacements')
    %         title(strcat('Map of interpolated displacements', ' @ ', dataName(1:4), ' mbar'));
    %         %axis([0 size_cell_image(2) 0 size_cell_image(1)]);
    %         pause(0.1)
    %     end
    %     figure
    %     %contourf(X/pix,Y/pix,d(i).dz_interp/pix,10)
    %     axis image
    
end

%% Calculate Q matrix. This is the matrix that relates tractions and displacements
% in Fourier space. The original derivation and expression for Q is in Xu
% et al. PNAS 107, 14964-14967 (2010). Note that this example code is for
% 2d traction force microscopy (where out-of-plane tractions are assumed to
% be negligible). If performing 3d TFM, ensure that the arguments to calcQ
% are modified appropriately. (see documentation for calcQ)

[nr,nc]=size(X);

fracpad=0;  %fraction of field of view to pad displacements on either side of current fov+extrapolated to get high k contributions to Q
nr2 = round((1+2*fracpad)*nr);
if mod(nr2,2)==0
    nr2=nr2+1;
end
Q = calcQ(h,h,EM,nu,nr2,dx,2); % Q matrix that interpolates between displacements and stresses at the substrate surface in Fourier space.

%% calculate filter for the displacement data (in Fourier space).
% This is effectively a low-pass exponential filter.
% No changes should be necessary if modifying code for 3d TFM.
min_feature_size=6;
qmax=nr2/(pi*min_feature_size);

% Get distance from of a grid point from the centre of the array
y=repmat((1:nr2)'-nr2/2,1,nr2);
x=y';
q=sqrt(x.^2+y.^2);

% Make the filter
qmsk=exp(-(q./qmax).^2);
%imagesc(qmsk)
qmsk=ifftshift(qmsk);

%% Calculate stresses from displacements
% If modifying code to perform 3d TFM, a corresponding z calculation needs
% to be added at each step.

% Make 1d Hann windows
[szr,szc]=size(d(1).dx_interp);
w_c=0.5*(1-cos(2*pi*(0:szc-1)/(szc-1)));
w_r=0.5*(1-cos(2*pi*(0:szr-1)/(szr-1)));

% Mesh Hann windows together to form 2d Hann window
[wnx,wny]=meshgrid(w_c,w_r);
wn=wnx.*wny;

%imagesc(wn)
%%

% Pad the window
padwidth=(nr2-nr)/2;
padheight=(nr2-nr)/2;
[sz1,sz2]=size(wn);
wn=[zeros(sz1+2*padheight,padwidth) [zeros(padheight,sz2);wn;zeros(padheight,sz2)] zeros(sz1+2*padheight,padwidth)];
% If you have the Image Processing Toolbox, this is equivalent to
% wn=padarray(wn,[(nr2-nr)/2,(nr2-nr)/2]);

for i = 1:length(d)
    % Get rid of NaN's in the interpolated displacement data
    d(i).dx_interp=extrapdisp(d(i).dx_interp);
    d(i).dy_interp=extrapdisp(d(i).dy_interp);
    %d(i).dz_interp=extrapdisp(d(i).dz_interp);
    % Pad and filter x displacements then multiply by the Hann window function
    [sz1,sz2]=size(d(i).dx_interp);
    utmp(1).u=d(i).dx_interp;%[zeros(sz1+2*padheight,padwidth) [zeros(padheight,sz2);d(i).dx_interp;zeros(padheight,sz2)] zeros(sz1+2*padheight,padwidth)];
    %If you have the Image Processing Toolbox, this line is equivalent to
    %utmp(1).u=padarray(d(i).dx_interp,[(nr2-nr)/2,(nr2-nr)/2]);
    utmp(1).u=real(ifft2(qmsk.*fft2(utmp(1).u)));
%     utmp(1).u %=utmp(1).u.*wn;
    
    % Pad and filter y displacements then multiply by the Hann window function
    utmp(2).u=d(i).dy_interp;%[zeros(sz1+2*padheight,padwidth) [zeros(padheight,sz2);d(i).dy_interp;zeros(padheight,sz2)] zeros(sz1+2*padheight,padwidth)];
    %If you have the Image Processing Toolbox, this line is equivalent to
    %utmp(2).u=padarray(d(i).dy_interp,[(nr2-nr)/2,(nr2-nr)/2]);
    utmp(2).u=real(ifft2(qmsk.*fft2(utmp(2).u)));
%     utmp(2).u %=utmp(2).u.*wn;
    
    %     % Pad and filter z displacements then multiply by the Hann window function
    %     utmp(3).u=[zeros(sz1+2*padheight,padwidth) [zeros(padheight,sz2);d(i).dz_interp;zeros(padheight,sz2)] zeros(sz1+2*padheight,padwidth)];
    %     %If you have the Image Processing Toolbox, this line is equivalent to
    %     %utmp(2).u=padarray(d(i).dy_interp,[(nr2-nr)/2,(nr2-nr)/2]);
    %     utmp(3).u=real(ifft2(qmsk.*fft2(utmp(3).u)));
    %     utmp(3).u=utmp(3).u.*wn;
    
    stmp = disp2stress(utmp,Q);
    
    % Remove the padding
    d(i).stress_x=stmp(1).s((nr2-nr)/2+1:((nr2-nr)/2+nr),(nr2-nr)/2+1:((nr2-nr)/2+nr));
    d(i).stress_y=stmp(2).s((nr2-nr)/2+1:((nr2-nr)/2+nr),(nr2-nr)/2+1:((nr2-nr)/2+nr));
    %d(i).stress_z=stmp(3).s((nr2-nr)/2+1:((nr2-nr)/2+nr),(nr2-nr)/2+1:((nr2-nr)/2+nr));
    d(i).ux=utmp(1).u((nr2-nr)/2+1:((nr2-nr)/2+nr),(nr2-nr)/2+1:((nr2-nr)/2+nr));
    d(i).uy=utmp(2).u((nr2-nr)/2+1:((nr2-nr)/2+nr),(nr2-nr)/2+1:((nr2-nr)/2+nr));
    %d(i).uz=utmp(3).u((nr2-nr)/2+1:((nr2-nr)/2+nr),(nr2-nr)/2+1:((nr2-nr)/2+nr));
    
    % Calculate the strain energy density = u.sigma/2 - see Mertz et al. PRL
    % 108, 198101 (2012).
    d(i).sed = 1/2*d(i).stress_x.*d(i).ux + 1/2*d(i).stress_y.*d(i).uy;
    
    % Calculate the stress magnitude
    d(i).mag = sqrt(d(i).stress_x.^2 + d(i).stress_y.^2);
    
    
    % At each timestep, plot the x/y displacements and the x/y traction
    % stresses in one figure. This will be empty at the reference time
    % step.
    %     figure
    %     subplot(2,1,1);
    %     imagesc([d(i).ux,d(i).uy])%,d(i).uz]);axis image;
    %     title(['time = ',num2str(i),', displacements']);
    %     subplot(2,1,2);
    % 	imagesc([d(i).stress_x,d(i).stress_y])%,d(i).stress_z]);axis image;
    % 	title('stresses');
    %     pause(1)
end
%% EXTRA IMAGES

% figure
% imagesc(X(1,:),Y(:,1),d(i).sed)%,d(i).stress_z]);axis image;
% %axis([0 size_cell_image(2)*pix 0 size_cell_image(1)*pix]);
% hold all
% %figure
% %imagesc(X(1,:),Y(:,1),sqrt(((d(i).stress_x.*d(i).fov).^2)+((d(i).stress_y.*d(i).fov).^2)))
% %quiver(X/pix,Y/pix,d(2).dx_interp/pix,d(2).dy_interp/pix)
% %title('stress mag');

%% Plot up various useful quantities
%
% Convert cell outline data from pixels to metres
%cell_x_metres=cell_x*pix;
%cell_y_metres=cell_y*pix;

% Plot up traction stress magnitude at the surface of the substrate
% (sigma.sigma).

% figure
% for i = 2
%    imagesc(X(1,:),Y(:,1),sqrt(((d(i).stress_x.*d(i).fov).^2)+((d(i).stress_y.*d(i).fov).^2))); colorbar;
%    hold on
%    plot(cell_x_metres,cell_y_metres,'w','LineWidth',2)
%    axis([0 size_cell_image(2)*pix 0 size_cell_image(1)*pix]);
%    hold off
%     pause(.1)
%     xlabel('x [m]')
%     ylabel('y [m]')
%     title('Traction stress magnitude [Pa]')
% end

%% Plot up strain energy density
%
% figure
% for i = 2
%     imagesc(X(1,:),Y(:,1),d(i).sed), colorbar
%     hold on
%     %plot(cell_x_metres,cell_y_metres,'w','LineWidth',2)
%     axis([0 size_cell_image(2)*pix 0 size_cell_image(1)*pix]);
%     hold off
%     xlabel('x [m]')
%     ylabel('y [m]')
%     title('Strain energy density [J/m^2]')
%     pause(0.1)
% end
% hold on;
% % plots velocities on top
%
% % temp = load(PIV);
% % quiver(4*temp.x*pix, 4*temp.y*pix, temp.u, temp.v, 'r');


%%
% Make quiver plots to show stresses on top of original image.

i=2;
%positions
str(i).r(:,1) = X(d(i).fov==1);
str(i).r(:,2) = Y(d(i).fov==1);
%stresses
str(i).dr(:,1) = d(i).stress_x(d(i).fov==1);
str(i).dr(:,2) = d(i).stress_y(d(i).fov==1);
sc = 4; %sc is scale for arrows
%  im = imcrop(im, rect./pix);
%size(str(i).mag)
size(str(i).r(:,1));
%figure
%imagesc(str(i).r(:,1)/pix,str(i).r(:,2)/pix,str(i).mag); hold on
imagesc(X(1,:)./pix, Y(:,1)./pix, d(i).mag), colorbar;
hold on;
axis([0 size_cell_image(2) 0 size_cell_image(1)]);
%colormap gray
quiver(str(i).r(:,1)/pix, str(i).r(:,2)/pix, sc*str(i).dr(:,1), sc*str(i).dr(:,2), 0, 'b');
%quiver(str(i).r(:,1), str(i).r(:,2), sc*d(i).stress_x, sc*d(i).stress_y, 0, 'b');
%hold off
xlabel('x [m]');
ylabel('y [m]');


% title(strcat('Stress Magnitude [N/m^2]', ' @ ', dataName(1:4), ' mbar'));
title(strcat('Stress Magnitude [N/m^2] @ ', iteration, ' s'));

%quiver(4*temp.x, 4*temp.y, temp.u, temp.v, 0.5, 'r');

% %% Scatterplot of all the displacements
% figure;
% scatter(str(i).dr(:, 1), str(i).dr(:, 2));
