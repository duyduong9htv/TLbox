%this class handles calculating the TL to a boxed area containing fish or
%scatters of interest. 
%The following properties have to be specified: 
% - location of source 
% - depth of source 
% -coordinates defining the region 
% -
classdef TLbox < handle 
    properties 
        zs; %source depth 
        xs; 
        ys; %coordinates (UTM) of source 
        nelements; 
        sourcetype; %source type, 'xf4' or 'mod30' in lower case
        x1; 
        y1; 
        x2; 
        y2; %coordinates of the vertices of rectangle box region of interest in UTM
        depth_start; 
        depth_end; %start and end depth for TL averaging
        depth_above_floor; %average TL only within, say 30 m off the sea bottom. 
        angle_step; %radial angle increment, in degrees
        transect; %transect (each radial direction, of to class TL) 
        TL; 
        f0; %frequency 
        L; %array aperture 
        MC; %number of Monte Carlo sample 
        rangeWindow; %TL averaging range window; 
        corr_length; %correlation length, = range to update sound speed profile to simulate the effect of internal wave
        RAMpath; %path to RAM executable file
        depthLimAve; %averaging within a depth window if set to 1.
        depthOffBottomAve; %averaging within some depths off the bottom if set to 1.
        radialTLangleAve;
        boxTLangleAve;
    end %properties 
    
    methods 
        %class constructor 
        function tlb = TLbox(xs, ys, zs, sourcetype, x1, y1, x2, y2, depth_start, depth_end)
            tlb.zs = 65; 
            
                     
            tlb.transect = TL; %make sure TL class in in the search path 
            tlb.xs = xs; 
            tlb.ys = ys; 
            tlb.zs = zs; 
            tlb.x1 = x1; 
            tlb.y1 = y1; 
            tlb.x2 = x2; 
            tlb.y2 = y2; 
            tlb.sourcetype = lower(sourcetype); 
            tlb.depth_start = depth_start; 
            tlb_depth_end = depth_end; 
            tlb.RAMpath = '/Users/dtran/RAM4/';
            tlb.MC = 3; 
                      
            
            %check source/receiver configuration: 
            if strcmp(tlb.sourcetype, 'mod30') %if MOD-30 source array is used 
                fid = fopen('source_spacing.in', 'wt'); 
                fprintf(fid,'10\n0.8382\n0.8001\n0.8382\n0.8001\n0.8001\n0.8382\n0.8382\n0.8382\n0.8382\n'); 
                fclose(fid); 
            elseif strcmp(tlb.sourcetype, 'xf4') %if XF-4 source array is used 
                fid = fopen('source_spacing.in', 'wt'); 
                fprintf(fid, '7\n1.6256\n1.6256\n1.6256\n1.6256\n1.6256\n1.6256\n'); 
                fclose(fid); 
            elseif strcmp(tlb.sourcetype, 'rcv')
                disp('Receiver configuration selected. Will use point source in RAM model!'); 
            else
                disp('Please select correct source/receiver type!'); 
                return;
            end
        end %end class constructor TLbox
        

        %select frequency and correct array aperture 
        function tlb = selectFrequency(tlb, f0) 
            tlb.f0 = f0; 
            tlb.transect.frequency = f0; 
            if f0 < 950 
                tlb.L = 94.5;
            else
                tlb.L = 47.25; 
            end
            
        end
        
        
        %calculate TL to a boxed area: 
        function tlb = getTLbox(tlb)
            
            initializeTLboxCalculation;             

            for ii = 1:size(radial_end,1) %loops through all the directions 
                disp(['Radial ' int2str(ii) ' out of ' num2str(size(radial_end,1))])
                initializeTransect; 
                
                for jj = 1:tlb.MC %looping through the number of Monte Carlo realizations
                    tlb.transect.randomSSP(); 
                    
                    if strcmpi(tlb.sourcetype, 'rcv')
                        tlb.transect.calculateGreen(); 
                    else
                        tlb.transect.calculateGreenSourceArray(); 
                    end 
                    
                    if jj==1
                        cf_mean_rangedep=abs(tlb.transect.gGrid).^2;
                        fine_depth=(dz:dz:dz*size(tlb.transect.gGrid,1));
                    else
                        cf_mean_rangedep=cf_mean_rangedep+abs(tlb.transect.gGrid).^2;
                    end
                        disp(['fraction done = ' num2str(jj/tlb.MC)]);
                end
                
                cf_mean_rangedep=cf_mean_rangedep/tlb.MC; %averaging

                %depth-averaging
                msq_fullda = zeros(nr,1);
                msq_eda = zeros(nr,1);
                water_depth = tlb.transect.bathy(1:nr,2);
                if tlb.depthLimAve == 1
                    msq_eda = aveWithinDepths(cf_mean_rangedep, water_depth, dz, tlb.depth_start, tlb.depth_end); %eda: effective depth averaged
                elseif tlb.depthOffBottomAve == 1
                    msq_eda = aveOffBottom(cf_mean_rangedep, water_depth, tlb.depth_above_floor, dz); 
                end
                
                %range-averaging
              
                range_ave_win=scale_win*log10(range_s);
                range_inc = dr;
                for rr=1:length(range_s)
                    if range_ave_win(rr)/2 > range_s(rr)     %situation where window would create negative index
                         ave_vec = (1:round((range_s(rr)+range_ave_win(rr)/2)/range_inc));
                    elseif (range_s(rr)+range_ave_win(rr)/2 > range_s(length(range_s)))  %situation where window would exceed bounds
                         ave_vec = (ceil((range_s(rr)-range_ave_win(rr)/2)/range_inc):1:length(range_s));
                    else
                         ave_vec = (ceil((range_s(rr)-range_ave_win(rr)/2)/range_inc):1:round((range_s(rr)+range_ave_win(rr)/2)/range_inc));
                    end
                    %average using hanning window in range:
                    han_win=hanning(length(ave_vec));
                    range_ave_han_eda(rr) = sum(han_win'.*msq_eda(ave_vec,1)')/sum(han_win);
                end

                TL_rda = 10*log10(range_ave_han_eda);
                datall(:,ii)=TL_rda';
                clear msq_fullda msq_eda TL_rda range_ave_han_eda range_ave_han_full cf_mean_rangedep;
                
            end %end looping in all direction 
     
            [radialTLangleAve, boxTLangleAve] = angleAverageTL(datall, az_ave_num, tlb.xs, tlb.ys, tlb.x1, tlb.x2, tlb.y1, tlb.y2, theta_vec, range_inc);
            tlb.radialTLangleAve = radialTLangleAve; 
            tlb.boxTLangleAve = boxTLangleAve; 
            
        end %end getTLbox function 
 
    end %methods 
end %classdef 
        
        