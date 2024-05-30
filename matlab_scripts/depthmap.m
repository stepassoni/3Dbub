function [depth_pos depth_neg depth_tot]=depthmap(BW, amplitude, pitch, phi, px_per_mm, deltaY, plot)

% Inputs:
% BW=frame to be processed (needed to calculate grid resolution
% amplitude=corrugation depth (2.5mm) 
% pitch (8.9mm), 
% angle= corrugation angle from flow direction in deg (-63 deg)
% plot= 0 or 1 depending wheter you want to plot the surface or not

%deltaY;
deltaX=0;

%[n_rows, n_cols] = size(BW);
n_rows=480;
n_cols=1280;
x_res=n_cols;
y_res=n_rows;

x_res_target=size(BW,2);
y_res_target=size(BW,1);

deltaX2=(x_res-x_res_target)/2;
deltaY2=(y_res-y_res_target)/2;

% Define the parameters
amplitude = amplitude/2; % mm
%pitch = 8.9; % mm
angle = phi; % degrees
angle2 = angle-180; % degrees

% Define the x and y ranges
Wd=n_cols/px_per_mm/2*cosd(angle)/pitch;

width=pitch/cosd(angle)*Wd/2;
height=y_res/px_per_mm;

x_range = linspace(0, width, x_res/4); % mm
y_range = linspace(0, height, y_res); % mm
[X, Y] = meshgrid(x_range, y_range);
Z = amplitude*sin(2*pi/pitch*(cosd(angle)*(X+deltaX) + sind(angle)*Y+deltaY))-amplitude;

x_range = linspace(width, 2*width, x_res/4); % mm
y_range = linspace(0, height, y_res); % mm
[X1, Y] = meshgrid(x_range, y_range);
Z1 = flip(Z,2);

x_range = linspace(0, width, x_res/4); % mm
y_range = linspace(0, height, y_res); % mm
[X, Y] = meshgrid(x_range, y_range);
Z2 = amplitude*sin(2*pi/pitch*(cosd(angle2)*(X+deltaX) + sind(-angle2)*Y+deltaY))+amplitude;

x_range = linspace(width, 2*width, x_res/4); % mm
y_range = linspace(0, height, y_res); % mm
[X1, Y] = meshgrid(x_range, y_range);
Z3 = flip(Z2,2);

if plot==1
    %Plot the surface
    hold on
    surf(X, Y, Z);
    surf(X1, Y, Z1);
    surf(X, Y, Z2);
    surf(X1, Y, Z3);
    surf(X+2*width, Y, Z);
    surf(X1+2*width, Y, Z1);
    surf(X+2*width, Y, Z2);
    surf(X1+2*width, Y, Z3);
    daspect([1 1 1])
    colormap;
    shading interp
    xlabel('X-axis [mm]');
    ylabel('Y-axis [mm]');
    zlabel('Z-axis [mm]');
    title('Corrugated Channel Depth');
    set(gca,'FontSize',14);
    colorbar()
end

depth=zeros(n_rows,n_cols);
depth_pos=zeros(n_rows,n_cols);
depth_neg=zeros(n_rows,n_cols);
for i=1:x_res/4
    depth(:,i)=Z(:,i)-Z2(:,i);
end    
for i=x_res/4+1:x_res/2
    depth(:,i)=Z1(:,i-x_res/4)-Z3(:,i-x_res/4);
end  
depth(:,x_res/2+1:end)=depth(:,1:x_res/2);

depth_pos(:,1:x_res/4)=Z2;
depth_pos(:,x_res/4+1:x_res/2)=Z3;
depth_pos(:,x_res/2+1:end)=flip(depth_pos(:,1:x_res/2),2);

depth_neg(:,1:x_res/4)=Z;
depth_neg(:,x_res/4+1:x_res/2)=Z1;
depth_neg(:,x_res/2+1:end)=flip(depth_neg(:,1:x_res/2),2);

depth_tot=abs(depth_pos-depth_neg);

depth_pos=flip(depth_pos);
depth_neg=flip(depth_neg);
depth_tot=flip(depth_tot);

depth_pos(:,1:deltaX2-1)=[];
depth_pos(:,n_cols-deltaX2*2+1:end)=[];
depth_pos(1:deltaY2-1,:)=[];
depth_pos(n_rows-deltaY2*2+1:end,:)=[];

depth_neg(:,1:deltaX2-1)=[];
depth_neg(:,n_cols-deltaX2*2+1:end)=[];
depth_neg(1:deltaY2-1,:)=[];
depth_neg(n_rows-deltaY2*2+1:end,:)=[];

depth_tot(:,1:deltaX2-1)=[];
depth_tot(:,n_cols-deltaX2*2+1:end)=[];
depth_tot(1:deltaY2-1,:)=[];
depth_tot(n_rows-deltaY2*2+1:end,:)=[];


% mask=depth_tot<0.2;
% BW(mask)=255;
% imshow(BW)