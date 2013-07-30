function [radialTLangleAve, boxTLangleAve] = ...
    angleAverageTL(datall, az_ave_num, src_utmx, src_utmy, x1, x2, y1, y2, theta_vec, range_inc)


    weight=hanning(az_ave_num);
    rep_weight=repmat(weight', size(datall,1), 1);
    datall_alog=10.^(datall/10);
    for az=1:size(datall,2)-az_ave_num+1
        dat_angle_ave_alog(:,az)=sum(datall_alog(:,az:az+az_ave_num-1).*rep_weight,2)/sum(weight);
    end
    radialTLangleAve=10*log10(dat_angle_ave_alog);
    theta_vec_ave=theta_vec((az_ave_num+1)/2:end-(az_ave_num-1)/2);
    dat_angle_ave = radialTLangleAve; 

    %% mapping radial to rectangular coordinate:
    grid_inc=30; 
    [X Y] = meshgrid(x1:grid_inc:x2, y2:-grid_inc:y1);
    R = sqrt((X-src_utmx).^2 +(Y-src_utmy).^2);
    THETA = atan2(Y-src_utmy,X-src_utmx);
    boxTLangleAve = interp2(theta_vec_ave,[0:length(dat_angle_ave)-1]*range_inc, dat_angle_ave, THETA, R);

end 