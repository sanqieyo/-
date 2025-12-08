function [x_length, y_length] = ProcessSeaArea( Param )
%ProcessSeaArea 根据视场范围自动生成海面大小
% 
y_list = Param.horizon_scale;           % 横向扫描序列
z_list = Param.vertical_scale;          % 纵向扫描序列
observe_dist = Param.observe_dist*1000;     %观测距离
rotate_factor = Param.rotate_factor;        %观测方向旋转矩阵
observe_position = Param.observe_position;      %观测点在目标坐标系下位置
x_0 = observe_position(1); 
y_0 = observe_position(2); 
z_0 = observe_position(3); 

y_loc_1 = observe_dist*tand(y_list(1));                   % 像素点在成像平面上的y坐标
z_loc_1 = tand(z_list(1))*observe_dist/cosd(y_list(1));      % 像素点在成像平面上的z坐标
scan_vector_1 = [-observe_dist, y_loc_1 , z_loc_1];        % 旋转后的坐标系中扫描射线的方向矢量
inc_vector_1 = scan_vector_1/rotate_factor;       % 扫描射线在原始坐标系的方向矢量

x_loc_1 = x_0 - inc_vector_1(1) / inc_vector_1(3) * z_0; 
y_loc_1 = y_0 - inc_vector_1(2) / inc_vector_1(3) * z_0;

m = length(y_list);
n = length(z_list);

y_loc_m = observe_dist*tand(y_list(m));                   % 像素点在成像平面上的y坐标
z_loc_1 = tand(z_list(1))*observe_dist/cosd(y_list(m));      % 像素点在成像平面上的z坐标
scan_vector_2 = [-observe_dist, y_loc_m , z_loc_1];        % 旋转后的坐标系中扫描射线的方向矢量
inc_vector_2 = scan_vector_2/rotate_factor;       % 扫描射线在原始坐标系的方向矢量
x_loc_m = x_0 - inc_vector_2(1) / inc_vector_2(3) * z_0; 

y_loc_1 = observe_dist*tand(y_list(1));                   % 像素点在成像平面上的y坐标
z_loc_n = tand(z_list(n))*observe_dist/cosd(y_list(1));      % 像素点在成像平面上的z坐标
scan_vector_3 = [-observe_dist, y_loc_1 , z_loc_n];        % 旋转后的坐标系中扫描射线的方向矢量
inc_vector_3 = scan_vector_3/rotate_factor;       % 扫描射线在原始坐标系的方向矢量
y_loc_n = y_0 - inc_vector_3(2) / inc_vector_3(3) * z_0; 

x_length = abs(round(2 * (x_loc_m - x_loc_1)));
if(mod(x_length, 2) == 1)
    x_length = x_length + 1;
end

y_length = x_length;

end

