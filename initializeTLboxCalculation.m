%parameters for running the TL problem 
load ssps_basin;
load ssps_bank;
c = 1500; 
lambda = c/tlb.f0; 
dr = 10*lambda; 
dz = lambda/10;             
Beta = 0.5*1.3*lambda/tlb.L; %broadside resolution
az_ave_num = 5; %number of transects to average in angle (use odd number)
max_depth=400;   % maximum depth allowed
corr_length = tlb.corr_length; 
% get range and angle limits 

[th1, r1] = cart2pol(tlb.x1 - tlb.xs, tlb.y1 - tlb.ys); 
[th2, r2] = cart2pol(tlb.x2 - tlb.xs, tlb.y1 - tlb.ys);
[th3, r3] = cart2pol(tlb.x1 - tlb.xs, tlb.y2 - tlb.ys);
[th4, r4] = cart2pol(tlb.x2 - tlb.xs, tlb.y2 - tlb.ys);

th = [th1 th2 th3 th4]; 
R = max([r1 r2 r3 r4]); 
th_min = min([th1 th2 th3 th4]); 
th_max = max([th2 th2 th3 th4]); 


az_ave_num=5;       %number of transects to average in angle (use odd number)
angle_overlap=Beta*(az_ave_num+1)/2;       %for azimuthal averaging
angle_dif=th_max-th_min;
if angle_dif>pi    %situation where box crosses negative y axis 
    for aa=1:4
        if th(aa)<0
            th(aa)=th(aa)+2*pi;
        end
    end
end
theta_vec=(th_min-angle_overlap):Beta:(th_max+angle_overlap);

%% range and depth

scale_win=1200;
range_overlap=scale_win*log10(R)/2;
range_max_RAM=R+range_overlap;

range_s =(dr:dr:range_max_RAM)';
nr = length(range_s);
range_ave_win=scale_win*log10(range_s);   %window for range averaging
rough_range=(0:corr_length:range_max_RAM);          %range sampling for ssps
radial_end=zeros(length(theta_vec),1);  
radial_end(:,1) = tlb.xs + range_max_RAM.*cos(transpose(theta_vec));  
radial_end(:,2) = tlb.ys + range_max_RAM.*sin(transpose(theta_vec));