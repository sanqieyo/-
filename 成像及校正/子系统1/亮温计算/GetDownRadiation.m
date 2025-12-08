function TB = GetDownRadiation(polar_vector,elevation_angle,Param)
% 此函数用于目标面元反射地面时的背景亮温
% 输入：极化向量polar_vector，入射向量俯仰角elevation_angle（俯为负），大气参数拟合系数targetdown_fitting_coeff，地面温度ground_temp
% 输出：入射方向上的环境辐射亮温TB

BRITEMP = polyval(Param.RadTemp.target_down_fitcoe(1,:),elevation_angle);           % 对应俯仰角下的目标与地面间的大气辐射
ATT = polyval(Param.RadTemp.target_down_fitcoe(2,:),elevation_angle);               % 对应俯仰角下的目标与地面间的大气衰减
BT_re = polyval(Param.RadTemp.target_down_fitcoe(3,:),elevation_angle);             % 对应俯仰角下地面反射方向上的天空亮温
[eh_ground,ev_ground] = RefGround(90+elevation_angle,Param.Ambient);            % 计算地面发射率分量
inc_vector = [-sind(90-elevation_angle),0,-cosd(90+elevation_angle)];          % 入射向量
normal_vector = [0,0,1];                        % 地面法向量
[e_ground,~] = LinearPolarization(polar_vector,inc_vector,normal_vector,eh_ground,ev_ground);     % 计地面反射率

ground_temp = Param.Ambient.ground_temp;
BT_ground = ground_temp*e_ground+(1-e_ground)*BT_re;        % 入射方向上的地面亮温
TB = BRITEMP +BT_ground*exp(-0.23*ATT);         % 考虑大气衰减和大气辐射
end