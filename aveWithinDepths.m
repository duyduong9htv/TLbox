function TLave = aveWithinDepths(cf_mean_rangedep, water_depth, dz, start_depth, end_depth)
% function TLave = aveWithinDepth(cf_mean_rangedep, water_depth, dz, start_depth, end_depth)


TLave = zeros(size(cf_mean_rangedep, 2)); 
for k = 1:size(cf_mean_rangedep, 2)
	start_ind = max(1, round(start_depth/dz)); 
	end_ind = min(round(water_depth(k)/dz), round(end_depth/dz)); 
	TLave(k) = mean(cf_mean_rangedep(start_ind:end_ind, k)); 
end
