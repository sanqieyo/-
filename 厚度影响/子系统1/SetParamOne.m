function SetParamOne( tt )
%SETPARAMONE 此处显示有关此函数的摘要
%   此处显示详细说明
global SimuParam;

%% 观测参数

Observe_matrix = [10 -68 90
                10 -70 90];
Move_matrix = [0 0 0
                0 0 0];

SimuParam.Observe.observe_dist = Observe_matrix(tt, 1);            % 观测距离（辐射计到目标），km
SimuParam.Observe.elevation_angle = Observe_matrix(tt, 2);         % 观测俯仰角（辐射计看目标），仰为正，角度制
SimuParam.Observe.azimuth_angle = Observe_matrix(tt, 3);           % 观测方位角（辐射计看目标，与x轴夹角），若y坐标为正则为正，角度制


%% 参数预处理
SimuParam.calculation_mode = 1;                           % 选择计算模式，0为常规计算，1为并行计算
[pixel_space,horizon_scale,vertical_scale] = ProcessScanParam(SimuParam.Observe);             % 扫描参数预处理
[rotate_factor,observe_position,radio_altitude] = ProcessCoordinate(SimuParam.Observe);            % 观测参数预处理
SimuParam.Observe.observe_position = observe_position;              % 观测点坐标
SimuParam.Observe.observe_position_1 = observe_position;

SimuParam.Observe.rotate_factor = rotate_factor;                    % 坐标旋转矩阵
SimuParam.Observe.radio_altitude = radio_altitude;                  % 辐射计高度
SimuParam.Observe.pixel_space = pixel_space;                        % 亮温分布成像平面间隔
SimuParam.Observe.horizon_scale = horizon_scale;                    % 亮温分布成像平面横向序列
SimuParam.Observe.vertical_scale = vertical_scale;                  % 亮温分布成像平面纵向序列

[x_length, y_length] = ProcessSeaArea(SimuParam.Observe);         %根据视场范围自动生成海面大小
SimuParam.Target.x_length = x_length;              %海面横向尺寸，单位米;
SimuParam.Target.y_length = y_length;              %海面纵向尺寸，单位米;

SimuParam.Target.sea.target_point = GetSeaPoint(SimuParam.Target,tt / 10);         % 目标节点坐标（x，y，z）
SimuParam.Target.sea.target_rank = GetSeaRank(SimuParam.Target);     % 目标面元（顶点编号）
SimuParam.Target.sea.target_temp = GetTemp(SimuParam.Target);

[sea_rank_point,sea_normal_vector,sea_target_point,shift_x,shift_y,shift_z] = ProcessSeaElement(SimuParam.Target,SimuParam.Observe, Observe_matrix, Move_matrix);                      % 面元参数预处理
[boat_rank_point,boat_normal_vector,boat_target_point,~,~,~] = ProcessBoatElement(SimuParam.Target,SimuParam.Observe, Observe_matrix, Move_matrix);                      % 面元参数预处理
SimuParam.Target.shift_x = shift_x;                                 % 目标坐标x的位移
SimuParam.Target.shift_y = shift_y;                                 % 目标坐标y的位移
SimuParam.Target.shift_z = shift_z;                                 % 目标坐标z的位移
SimuParam.Target.sea.rank_point = sea_rank_point;                           % 按面元节点排列的（x,y,z）坐标
SimuParam.Target.sea.target_point = sea_target_point;
SimuParam.Target.sea.normal_vector = sea_normal_vector;                     % 面元法向量
SimuParam.Target.boat.normal_vector = boat_normal_vector;         % 目标节点坐标（x，y，z）
SimuParam.Target.boat.target_point = boat_target_point;
SimuParam.Target.boat.rank_point = boat_rank_point;     % 目标面元（顶点编号）\ 


polar_vector = ProcessPolarization(SimuParam.Observe);                %极化参数预处理
SimuParam.Observe.polar_vector = polar_vector;                      % 极化向量


[radio_atmosparam,target_atmosparam,targetdown_angle,dist_radiation,dist_attenuation] =AtmosInitialize(SimuParam);             % 大气参数预处理


SimuParam.RadTemp.target_atmosparam = target_atmosparam;            % 目标处背景亮温
SimuParam.RadTemp.radio_atmosparam = radio_atmosparam;              % 辐射计处背景亮温
SimuParam.RadTemp.targetdown_angle = targetdown_angle;              % 目标处向下观测时反射地面临界角
SimuParam.RadTemp.dist_radiation = dist_radiation;                  % 目标与辐射计间的大气辐射
SimuParam.RadTemp.dist_attenuation = dist_attenuation;              % 目标与辐射计间的大气衰减

end

