function [ boat_target_point, sea_target_point] = ProcessMove( Param, a, b )
%PROCESSMOVE 此处显示有关此函数的摘要
%   此函数用于解决辐射计在成像过程中的移动轨迹问题
%   将原有的position中的原点移动到新的成像中心
boat_target_point = Param.Target.boat.target_point;
sea_target_point = Param.Target.sea.target_point;
observe_position = Param.Observe.observe_position;


theta = -atan(observe_position(1,2) / observe_position(1,1));
rotate_factor = [cos(theta) sin(theta) 0
                 -sin(theta) cos(theta) 0
                 0 0 1];
observe_position = observe_position * rotate_factor - b;
observe_position = observe_position / rotate_factor;

rotate_factor_1 = Param.Observe.rotate_factor;
inc_vector = [a(1,1)*1000 0 0];
inc_vector = inc_vector / rotate_factor_1;
x_0 = observe_position(1);
y_0 = observe_position(2);
z_0 = observe_position(3);
x_loc = x_0 - inc_vector(1) / inc_vector(3) * z_0; 
y_loc = y_0 - inc_vector(2) / inc_vector(3) * z_0;

boat_target_point(:,1) = boat_target_point(:,1) - x_loc;
boat_target_point(:,2) = boat_target_point(:,2) - y_loc;

sea_target_point(:,1) = sea_target_point(:,1) - x_loc;
sea_target_point(:,2) = sea_target_point(:,2) - y_loc;

% observe_position = (observe_position * rotate_factor + b) / rotate_factor;




end

