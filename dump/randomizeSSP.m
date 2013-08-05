shidx = find(bathy_data(:,2)<=80);
if ~isempty(shidx)
    basin = floor(bathy_data(shidx(1),1)/corr_length);
    rand_ssps = ceil((size(ssps_basin,2)-1).*rand(basin,1))+1;
    rand_ssps = [rand_ssps;ceil((size(ssps_bank,2)-1).*rand(length(rough_range)-basin,1))+1];
    for ss = 1:basin
        rough_svp(:,ss) = ssps_basin(:,rand_ssps(ss));
    end
    for ss = basin+1:length(rough_range)
        rough_svp(:,ss) = ssps_bank(:,rand_ssps(ss));
    end
else
    rand_ssps=ceil((size(ssps_basin,2)-1).*rand(length(rough_range),1))+1; 
    for ss=1:length(rough_range)
        rough_svp(:,ss)=ssps_basin(:,rand_ssps(ss));
    end
end
