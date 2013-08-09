% This code demonstrates how to use the TLbox class to calculate the
% transmission loss from either a source array to a rectangular region or
% from a rectangular region to a receiver. 


%% initialize 

% tlb = TLbox(xs, ys, zs, sourcetype, x1, y1, x2, y2, depth_start, depth_end)
tlb = TLbox(608202, 4673791, 65, ...
            'rcv', ...
            608202+10e3, 4673791-20e3, ...
            608202+30e3, 4673791-10e3, ...
            100, 200); 
        
%options for source type includes: 'xf4', 'mod30', or 'rcv'. Choosing 'xf4'
%or 'mod30' will enable source array mode. Choosing 'rcv' will only use
%point source TL calculation. 

tlb.MC = 1; 
tlb.depthOffBottomAve = 1; %set the mode to averaging within some water column off of the bottom 
tlb.depth_above_floor = 50; %water column of 50 m from the bottom 


%other options: 
% tlb.depthLimAve = 1; %set the averaging mode to say from 100 m to 150 m.
% This depends on where the fish are hypothesized/found to be. If
% tlb.depthLimAve is chosen, remember to set tlb.depth_start and
% tlb.depth_end. If tlb.depth_end goes beyond the sea bottom, the program
% will limit the averaging window only up to the bottom interface. 


tlb.selectFrequency(415); %set the frequency to 415 Hz. The aperture L will be set accordingly 
tlb.corr_length = 500; %range to update the sound speed profile 
tlb.setRAMpath('/Users/dtran/Research/TransmissionLoss/TL/fortran/MacOS/'); %path where the executable RAM files are located. 

tlb.getTLbox; %calculates TL to the boxed area, loop through all radial directions 




%results check: 
figure; imagesc([tlb.x1 tlb.x2], [tlb.y2 tlb.y1], tlb.boxTLangleAve); caxis([-80 -70]); axis xy
hold on; 



%% NOTES 


% The default number of Monte Carlo runs is set to 3. The user is free
% to change to include more randomization and averaging to smooth the TL. 


% The program uses the bathymetry from the Gulf of Maine that is stored in
% the file bathymetry.arr. If other environment was used, please make sure
% the transect bathymetries are correct. Replace the function
% getTransectUTM in the TL class by your own transect creating function 





