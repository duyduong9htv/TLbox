function TLave = aveOffBottom(cf_mean_rangedep, water_depth, depth_above_floor, dz)
% function TLave = aveOffBottom(cf_mean_rangedep, water_depth, depth_above_floor, dz)


TLave = zeros(size(cf_mean_rangedep, 2)); 
for k = 1:size(cf_mean_rangedep, 2)
	start_ind = max(1, round(depth_above_floor/dz)); 
	end_ind = round(water_depth(k)/dz); 
	TLave(k) = mean(cf_mean_rangedep(start_ind:end_ind, k)); 
end


	