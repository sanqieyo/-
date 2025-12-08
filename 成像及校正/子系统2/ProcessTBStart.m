function [ TB_large, TB_original,TB_error] = ProcessTBStart( Param, TB_original)
% 用正负5度小视场输入图得到一张120度的图

% 输入：     TB_original   子系统1输出的亮温矩阵

% 输出：     TB_large      大的背景亮温矩阵
%            TB_original   原始亮温图-均匀背景

fov_angle = Param.Scene.fov_angle;
FOP_theta = Param.AntPattern.FOP_theta;  
fov_large = FOP_theta * 2;
theta_original = linspace(-fov_angle/2 + 90, fov_angle/2 + 90,size(TB_original,1));
phi_original = linspace(-fov_angle/2 + 90, fov_angle/2 + 90,size(TB_original,1));
[Theta_original,Phi_original] = ndgrid(theta_original,phi_original);

% theta_large = linspace(-FOP_theta + 90, FOP_theta + 90,5*size(TB_original,1));
% phi_large = linspace(-FOP_theta + 90, FOP_theta + 90,5*size(TB_original,1));
theta_large = linspace(-FOP_theta + 90, FOP_theta + 90,size(TB_original,1));
phi_large = linspace(-FOP_theta + 90, FOP_theta + 90,size(TB_original,1));
[Theta_large,Phi_large] = ndgrid(theta_large,phi_large);

TB_avg = (sum(TB_original(1,:)) + sum(TB_original(size(TB_original,1),:)) + sum(TB_original(:,1)) + sum(TB_original(:,size(TB_original,1))))/size(TB_original,1)/4;

TB_large = TB_avg*ones(length(theta_large), length(phi_large)); % 均匀背景值 = 原始亮温最外一圈数据均值
TB_original = TB_original - TB_avg; % 原始图-均匀背景

% TB_tar_scene  = griddata(theta_large,phi_large,scene_ori,theta_original, phi_original); 
% TB_original = TB_original - TB_tar_scene;
% TB_large = scene_ori;

% TB_large = griddata(Theta_original,Phi_original,TB_original,Theta_large, Phi_large); 
% TB_value = TB_large(~isnan(TB_large));
% % average_TB = sum(TB_value)/length(TB_value);
% TB_large(~isnan(TB_large)) = 0;
% TB_large(isnan(TB_large)) = 230;

theta_error = linspace(-fov_angle/2 + 90 - 0.5, fov_angle/2 + 90 + 0.5,51);
phi_error = linspace(-fov_angle/2 + 90 - 0.5, fov_angle/2 + 90 + 0.5,51);
[Theta_error,Phi_error] = ndgrid(theta_error,phi_error);
TB_error  = griddata(Theta_original,Phi_original,TB_original,Theta_error, Phi_error); 

end