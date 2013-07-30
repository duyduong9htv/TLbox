% This file calculates the transmission loss from a source to a rectangular
% area specified by the user (usually a box of interest with significant
% density, in order to get the final scatter strength map) 

%INPUTS: 
% -point source location (will be vertical array source if needed)
% -boxed region of interest (ROI) 
% -location of pixels within ROI 

%METHOD: 
% -calculate the 2 radial directions Phi_1, Phi_2 that bound the ROI
% 
% -get the bathymetry and calculate transmission losses along all
% directions between Phi_1 and Phi_2 at an angle step size equal to the
% array broadside resolution (best angle resolution affordable with the
% array). This step size can be changed, if needed. 
% 
% -convert polar coordinates to cartesian coordinates, interpolate for all
% values within the ROI. 
%
% -Averaging: the TL should be range-averaged and averaged across the angle
% OR: the TL can just be calculated (interpolated) for all points within
% the ROI, then do a 2-D smoothing. 
% 
% -this is a rewritten version of Mark's previous TL calculation code. 
% Improvements are: 
%   -remove dependencies (have to load a bunch of dependence files *.in,
%   INFILErcv*, *.m etc) by making full use of the database structure query
%   -make use of the TL class to calculate one way TL in each transect, get
%   the functions better organized and reduce the length of the code. 


cd /Volumes/Neptune2/Duong/WHALES/whale_localization_data/Tracks_data/track570_4/ARR

filename = 'fora2006jd276t023615'; 
pings = PingQuery('track', '570_4'); %ping index: 126 for the corresponding filename 


%% Initial parameters 

f = 950; %frequency 
c = 1500; %speed of sound (constant) 
lambda = c/f; 
L = 47.5; %array aperture 

up_depth = 80; % upper bound of fish school depth 
low_depth = 100; %lower depth of fish school 

angle_step = 1.5*lambda/L; %broad side angular resolution 
dr = 30; %range increment 
dz = 0.2; %depth increment 

load fish_box; %rectangle ROI 
srcPos = pings(10).srcUTM; %location of source for PE runs
MC = 3; %number of Monte Carlo simulations to run 


%% calculation

% r = [ddist(srcPos, [fish_box(1) fish_box(3)]); ...
%      ddist(srcPos, [fish_box(2) fish_box(3)]); ...
%      ddist(srcPos, [fish_box(1) fish_box(4)]); ...
%      ddist(srcPos, [fish_box(2) fish_box(4)])]; %range from point source (receiver) to all corners of ROI 
 

[th1, r1] = cart2pol(fish_box(1) - srcPos(1), fish_box(3) - srcPos(2)); 
[th2, r2] = cart2pol(fish_box(2) - srcPos(1), fish_box(3) - srcPos(2));
[th3, r3] = cart2pol(fish_box(1) - srcPos(1), fish_box(4) - srcPos(2));
[th4, r4] = cart2pol(fish_box(2) - srcPos(1), fish_box(4) - srcPos(2));

R = max([r1 r2 r3 r4]); 
th_min = min([th1 th2 th3 th4]); 
th_max = max([th2 th2 th3 th4]); 


r = dr:dr:R; 
theta = th_min:angle_step:th_max; 

radial_grid = zeros(length(r), length(theta)); %initialize the TL grid for all r and theta values

th_count = 0; 
tic 
for th = th_min:angle_step:th_max %sweep all the directions 
    th_count = th_count +1;
    disp((th - th_min)/(th_max - th_min)); 
    %set up the bathymetry 
    [x, y] = pol2cart(th, R); 
    x = x + srcPos(1); 
    y = y + srcPos(2); 
    rayTL = TL; 
    rayTL.frequency = f; 
    rayTL.maxRange = R + 1e3; 
    rayTL.zmax = 300; 
    rayTL.ranges = 0:500:rayTL.maxRange; 
    rayTL.zs = 65; %from source array 
    rayTL.x1 = srcPos(1); rayTL.y1 = srcPos(2); 
    rayTL.x2 = x; rayTL.y2 = y; 
    rayTL.dr = dr;  
    rayTL.getTransectUTM(); %get the bathymetry in this direction 
    
    %% MC simulations 
    lineTL = 0; 
    for k = 1:MC %loop through number of Monte Carlo simulations 
        disp(k); 
        rayTL.randomSSP; 
        rayTL.calculateGreen; 

        %depth-averaging of the TL over the fish school depth: 
        depth_inds = round(up_depth/rayTL.dz):1:round(low_depth/rayTL.dz); 
        lineTL = lineTL + mean(abs(rayTL.gGrid(depth_inds, :)).^2); 
    end
    
    lineTL = lineTL/MC; 
    lineTL = lineTL(1:length(r)); 
    
    radial_grid (:, th_count) = 10*log10(lineTL);
    
end

% interpolate for ROI TL 

xx = fish_box(1):30:fish_box(2); 
yy = fish_box(3):30:fish_box(4); 
TLbox = zeros(length(xx), length(yy)); 
radial_grid2 = radial_grid; 
radial_grid = smooth2a(radial_grid2, 5, 5); 

for kx = 1:length(xx)
    disp(kx/length(xx)); 
    for ky = 1:length(yy)
        [t1, r1] = cart2pol(xx(kx) - srcPos(1), yy(ky) - srcPos(2)); 
        ind_th = find(abs(theta(:) - t1) == min(abs(theta(:) - t1))); 
        ind_r = find(abs(r(:) - r1) == min(abs(r(:) - r1)));         
        TLbox(kx, ky) = radial_grid(ind_r, ind_th); 
    end
end
toc 
