function [ksai_cur, eta_cur, offset_cur, theta_nd, phi_nd] = GetOriginalLM( TB_cope,sc_fov_angle )
% 输入背景亮温矩阵和天线观测角，输出L、M，倾斜因子和theta,phi
% 输入：     TB_cope 背景矩阵
%            sc_fov_angle 天线观测角（如果是正负60，这里就是120）
% 输出：     ksai_cur 即l = sin(theta)*cos(phi)
%            eta_cur 即m = sin(theta)*sin(phi)
%            offset_cur 倾斜因子，即sqrt(1-l^2-m^2)

theta_range = linspace(-sc_fov_angle/2 + 90, sc_fov_angle/2 + 90,size(TB_cope,1));
phi_range = linspace(sc_fov_angle/2, -sc_fov_angle/2,size(TB_cope,1));
[theta_nd, phi_nd] = ndgrid(theta_range,phi_range);
ksai_cur = sind(theta_nd).*cosd(phi_nd);
eta_cur = sind(theta_nd).*sind(phi_nd);

offset_cur = real(sqrt(1 - ksai_cur.*ksai_cur - eta_cur.*eta_cur ));
z_zhengfu = -ones(size(theta_nd));                    % 全部置为-1
z_zhengfu(ceil(size(theta_nd,1) / 2),:) = 0;          % 251行是0
z_zhengfu(1 : ceil(size(theta_nd,1) / 2) - 1, :) = 1; % 前250行都是1，后250行都是-1
offset_cur = offset_cur.*z_zhengfu;

end



