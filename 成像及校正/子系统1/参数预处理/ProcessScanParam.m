function [pixel_space,horizon_scale,vertical_scale]=ProcessScanParam(Param)
% 本函数用于对亮温分布成像平面初始化
% 输入横向扫描次数horizontal_scan，纵向扫描次数vertical_scan，扫描角度步长scan_step，亮温分布图像精度scan_factor
% 输出成像平面像素点的距离间隔pixel_space、横向像素范围horizontal_scale、纵向像素范围vertical_scale、像素点坐标pixel_position

Ymax = Param.scan_step*Param.horizon_scan/2;         % 由视场范围得出成像平面区域
Zmax = Param.scan_step*Param.vertical_scan/2;        % 将角度范围转换为距离范围
pixel_space = Param.scan_step*Param.scan_factor;         % 由扫描角度步长换算为成像平面像素点间隔
Ymax = ceil(Ymax/pixel_space)*pixel_space;           % 将Ymax取整为pixel_space的整数倍，确保成像平面两端对称
Zmax = ceil(Zmax/pixel_space)*pixel_space;           % 将Zmax取整为pixel_space的整数倍，确保成像平面两端对称
horizon_scale = (-Ymax:pixel_space:Ymax);            % 将成像平面离散化
vertical_scale = (Zmax:-pixel_space:-Zmax);
end