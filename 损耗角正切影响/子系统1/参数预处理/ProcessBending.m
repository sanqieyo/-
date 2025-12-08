function [ sea_point ] = ProcessBending( pre_sea_point, R)
%PROCESSBENDING 利用地球半径对海面进行几何弯曲

sea_point = ones(length(pre_sea_point(:, 1)),3);
for ii = 1 : length(pre_sea_point(:, 1))
    coordinate = pre_sea_point(ii, :);
    x_0 = coordinate(1,1);
    y_0 = coordinate(1,2);
    z_0 = coordinate(1,3);
    
    lr = sqrt(x_0*x_0 + y_0*y_0);
    alpha = atand(x_0 / y_0);                               %绕z轴弯曲角度alpha
    beta = acosd((2 * R * R - lr * lr) / (2 * R * R));     %绕y轴弯曲角度beta
    
    rotate_1 = [cosd(alpha) -sind(alpha) 0
                sind(alpha) cosd(alpha) 0 
                0 0 1];
    rotate_2 = [cosd(beta) 0 sind(beta)
                0 1 0
                -sind(beta) 0 cosd(beta)];
    rotate = rotate_1 * rotate_2;
    
    sea_point(ii,:) = [lr, 0 , z_0] / rotate ;
end

