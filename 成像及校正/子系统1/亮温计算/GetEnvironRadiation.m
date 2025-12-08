function TB = GetEnvironRadiation(polar_vector,elevation_angle,Param,inc_vector)
% 函数名：GetEnvironRadiation
% 功能：当目标面元反射环境时，计算反射方向上的环境辐射
% 输入：入射方向上的极化向量polar_vector，反射向量俯仰角elevation_angle，多次反射次数trace_count，入射向量inc_vector
% 输出：反射方向上的环境辐射TB
index = ceil((elevation_angle-(-90)+0.00001)/0.1);              % 判断俯仰角所在区间
factor = ceil((elevation_angle-(-90)+0.00001)/0.1)-(elevation_angle-(-90)+0.00001)/0.1;                        % 加权因子
target_atmos = Param.RadTemp.target_atmosparam;                 % 目标处背景亮温
if elevation_angle < Param.RadTemp.targetdown_angle
    BRITEMP = target_atmos(1,index)*factor+target_atmos(1,index+1)*(1-factor);           % 对应俯仰角下的目标与地面间的大气辐射
    ATT = target_atmos(2,index)*factor+target_atmos(2,index+1)*(1-factor);               % 对应俯仰角下的目标与地面间的大气衰减
    BT_re = target_atmos(3,index)*factor+target_atmos(3,index+1)*(1-factor);             % 对应俯仰角下地面反射方向上的天空亮温
    [eh_ground,ev_ground] = RefGround(90+elevation_angle,Param.Ambient);            % 计算地面发射率分量
    normal_vector = [0,0,1];                        % 地面法向量
    e_ground = LinearPolarization(polar_vector,inc_vector,normal_vector,eh_ground,ev_ground);     % 计地面反射率
    ground_temp = Param.Ambient.ground_temp;
    BT_ground = ground_temp*e_ground+(1-e_ground)*BT_re;        % 入射方向上的地面亮温
    TB = BRITEMP +BT_ground*exp(-0.23*ATT);         % 考虑大气衰减和大气辐射
    
else
    TB = target_atmos(1,index)*factor+target_atmos(1,index+1)*(1-factor);        % 加权计算对应角度下的目标面元温度
end
end