function [rotate_factor,position,radio_altitude]=ProcessCoordinate(Param)
%本函数用于计算目标坐标系的旋转矩阵
%输入结构体SimuParam
%输出旋转矩阵rotate_factor，观测点在目标坐标系中的坐标position，目标离地面高度target_altitude

R = 6378.14;                % 地球半径，单位km
%% 不考虑地球半径

% % 计算旋转矩阵
% a=sind(-Param.elevation_angle);
% b=cosd(-Param.elevation_angle);
% c=sind(-Param.azimuth_angle);
% d=cosd(-Param.azimuth_angle);
% rotate_factor=[b*d (b^2+a^2)*c -a*d
%               -b*c d a*c
%              (a^3+a*b^2)*(1-d)+a*d 0 (b^3+b*a^2)*(1-d)+b*d]; %坐标转换矩阵
% elevation_angle = Param.elevation_angle+90;       % elevation_angle为观测矢量与目标z轴夹角，azimuth_angle为观测矢量与目标x轴夹角
% position = [Param.observe_dist*sind(elevation_angle)*cosd(Param.azimuth_angle),Param.observe_dist*...
%           sind(elevation_angle)*sind(Param.azimuth_angle),Param.observe_dist*cosd(elevation_angle)];     %计算观测点在目标坐标系中的位置
% position = position*1000;           % 将单位转换为m
% radio_altitude = Param.target_altitude+Param.observe_dist*cosd(elevation_angle);                                 %计算目标离地面高度

%% 考虑地球半径
% 考虑地球半径，求解辐射计高度
D = Param.observe_dist;                 % 观测距离，单位km
y = Param.target_altitude+R;            % 目标离地心距离，单位km
theta = Param.elevation_angle+90;       % 观测矢量与辐射计到球心方向矢量夹角，单位角度
x =  D*cosd(theta) + sqrt(y^2 - D^2 * sind(theta)^2);              % 辐射计距球心距离(杨辰修改过的)
radio_altitude = x-R;                   % 辐射计离地面高度

% 考虑地球半径，计算辐射计在目标坐标系中的坐标
elevation_angle = real(acosd((D^2+x^2-y^2)/(2*D*x)));                   % 忽略虚部
% elevation_angle为观测矢量与目标z轴夹角，azimuth_angle为观测矢量与目标x轴夹角
position = [Param.observe_dist*sind(elevation_angle)*cosd(Param.azimuth_angle),Param.observe_dist*...
          sind(elevation_angle)*sind(Param.azimuth_angle),Param.observe_dist*cosd(elevation_angle)];     %计算观测点在目标坐标系中的位置
position = position*1000;           % 将单位转换为m

% 计算旋转矩阵
elevation_angle = elevation_angle-90;
a=sind(-elevation_angle);
b=cosd(-elevation_angle);
c=sind(-Param.azimuth_angle);
d=cosd(-Param.azimuth_angle);
rotate_factor=[b*d (b^2+a^2)*c -a*d
              -b*c d a*c
             (a^3+a*b^2)*(1-d)+a*d 0 (b^3+b*a^2)*(1-d)+b*d]; %坐标转换矩阵
end