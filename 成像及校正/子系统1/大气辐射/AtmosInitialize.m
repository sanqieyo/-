function [radio_atmosparam,target_atmosparam,targetdown_angle,dist_radiation,dist_attenuation] = AtmosInitialize(Param)
% 本函数用于大气辐射计算的初始化
% 输入结构体SimuParam，目标高度Param
% 输出辐射计背景亮温拟合系数radio_fitting_coeff，目标背景亮温拟合系数target_fitting_coeff
% 以及辐射计与目标间的大气辐射Dist_Tsky，大气衰减因子Loss
Param = ProcessAtmosParam(Param);   % 建立大气分层模型

%% 计算背景亮温

[radio_atmosparam,target_atmosparam,targetdown_angle] = GetBackground(Param);      % 计算辐射计及目标背景亮温

%% 计算大气衰减

[BRITEMP,ATT,mode] = GetAtmosAttenuation(Param);    % 计算辐射计与目标间的大气衰减和大气辐射
if mode == 1
   dist_attenuation = exp(0.23*ATT);    % 存储衰减率
   dist_radiation = BRITEMP;            % 存储辐射计和目标间的大气亮温
else if mode ==0
        error('无法探测到目标！！');
    end
end
end
