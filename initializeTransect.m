%%initialize a transect TL problem 
tlb.transect.frequency = tlb.f0; 
tlb.transect.zmax = 400; 
tlb.transect.maxRange = range_max_RAM; 
tlb.transect.ssps_bank = ssps_bank; 
tlb.transect.ssps_basin = ssps_basin; 
tlb.transect.dr = dr; 
tlb.transect.dz = dz; 
tlb.transect.zs = tlb.zs; 
tlb.transect.zr = 105; 
tlb.transect.x1 = tlb.xs; 
tlb.transect.y1 = tlb.ys; 
tlb.transect.x2 = radial_end(ii, 1); 
tlb.transect.y2 = radial_end(ii, 2); 
%                 tlb.transect.RAMpath = tlb.RAMpath; 

                
tlb.transect.getTransectUTM();