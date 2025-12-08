function SetParam()
% 用于仿真程序的参数设置以及仿真参数的预处理
global SimuParam;

%% 观测参数

SimuParam.Observe.target_altitude = 0;          % 目标离地面高度，km
SimuParam.Observe.horizon_scan = 100;        % 辐射计横向扫描次数
SimuParam.Observe.vertical_scan = 100;          % 辐射计纵向扫描次数
% SimuParam.Observe.scan_step = 0.02;              % 辐射计扫描步长，角度制
SimuParam.Observe.scan_step = 0.1; 
SimuParam.Observe.scan_factor = 1;            % 亮温分布图像精度，参数范围0-1，越接近于0精度越高
SimuParam.Observe.polar_angle = 0;              % 辐射计极化角，0为水平极化，90为垂直极化，角度制
SimuParam.Observe.beam_width = 0.5;     % 辐射计天线波束宽度，角度制
SimuParam.Observe.center_freq = 94;             % 辐射计中心频率，GHz

%% 目标参数

SimuParam.Target.pitch_angle = 0;                       % 目标飞行姿态-俯仰角，角度制
SimuParam.Target.yaw_angle = 0;                         % 目标飞行姿态-偏航角，角度制
SimuParam.Target.roll_angle = 0;                        % 目标飞行姿态-翻滚角，角度制
SimuParam.Target.target_temp = 293;      % 目标物理温度，K
SimuParam.Target.material_type = 1;                      % 目标表面材料类型 ，0为金属，1为隐身涂层（单层模型），2为隐身涂层（双层模型）
SimuParam.Target.target_conductivity = 5*10^7;          % 目标金属电导率                                                         
SimuParam.Target.layer1_target_permittivity = 8.32 - 13.68i;        % 第1层（红外），目标介电常数，红外2.3-0.02i，雷达20.6-3.07i，雷达18.4-2.57i
SimuParam.Target.layer1_target_permeability = 1;       % 第1层，目标表面涂层磁导率，红外0.67-0.75i，雷达1.24-0.93i，雷达1.03-0.72i
SimuParam.Target.layer1_coat_thickness =1000;              % 第1层，目标表面涂层厚度（单位m），红外0.000033，雷达0.0005，雷达0.0004


SimuParam.Target.boat.target_rank = load('模型数据\deng2_rank_3.txt');     % 目标面元（顶点编号）
SimuParam.Target.boat.target_point = load('模型数据\deng2_point_3.txt');         % 目标节点坐标（x，y，z）
SimuParam.Target.boat.target_temp = load('模型数据\deng2_temp_3.txt');


SimuParam.Target.wind_angle = 0;             %海面上空风向角,单位 度;
SimuParam.Target.wind_speed = 10;            %海面上空风速 U_10,单位m/s;
SimuParam.Target.delta_x = 0.5;              %海面横向最小，单位米;
SimuParam.Target.delta_y = 0.5;              %海面横向尺寸，单位米;

 
% SimuParam.Target.target_rank = GetSeaRank;     % 目标面元（顶点编号）
% SimuParam.Target.target_point = GetSeaPoint;         % 目标节点坐标（x，y，z）
% SimuParam.Target.target_temp = GetTemp;
%% 环境参数
SimuParam.Ambient.atmos_permittivity = 8.8e-12;         % 自由空间介电常数
% SimuParam.Ambient.ground_permittivity = 6.1955 - 0.3386i;       % 地面（或海面）介电常数
SimuParam.Ambient.ground_permittivity = 8.32 - 13.68i;       % 地面（或海面）介电常数
SimuParam.Ambient.ground_temp = 295;            % 地面（或海面）物理温度
SimuParam.Ambient.atmos_type = 'SUMMER MIDLAT';                 % 可选项：'SUMMER MIDLAT'，'WINTER MIDLAT','SUMMER SUBART','WINTER SUBART','ANNUAL TROPIC','ITU-R REF STD'
SimuParam.Ambient.weather_type = 'FAIR';                        % 可选项：'FAIR'晴空大气，'CLOUDY'云，'MIST '雾，'RAINY '微雨
SimuParam.Ambient.earth_radius = 6378.14;                       % 地球半径，单位km


%% 系统参数设置
SimuParam.SystemInput.Tau = 0.1; % 积分时间
SimuParam.SystemInput.c = 3.0e8; % 光速
SimuParam.SystemInput.f0 = 94e9; % 中心频率
SimuParam.SystemInput.Ta = 265;  % 天线温度
SimuParam.SystemInput.Tr = 434;  % 接收机等效噪声温度
SimuParam.SystemInput.BandWidth = 2000.e6; % 带宽
SimuParam.SystemInput.gain = 23;           % 天线增益/dB，方向性系数

%% 天线方向图
SimuParam.AntPattern.pat_file_path = '子系统2\参数预处理\参数\94new.txt';
SimuParam.AntPattern.FOP_theta = 60;
SimuParam.AntPattern.pat_pix_count = 501;
SimuParam.AntPattern.polarization = [0 1 0]; % 天线接收电场方向

%% 阵列设置
SimuParam.Array.array_mat = '子系统2\参数预处理\参数\array.mat';
% SimuParam.Array.array_mat = '子系统2\参数预处理\参数\array_Y.mat';
SimuParam.Array.del_u_input = 0.021; %六边形
% SimuParam.Array.del_u_input = 0.02; % Y

%% 反演视场设置
SimuParam.InvFov.Scene_Xi_inv = sind(5)*2/sqrt(3);% 反演视场范围
SimuParam.InvFov.Scene_Eta_inv = sind(5)*2/sqrt(3);
% SimuParam.InvFov.Scene_Xi_inv = sind(7/3)*2/sqrt(3);% 反演视场范围
% SimuParam.InvFov.Scene_Eta_inv = sind(7/3)*2/sqrt(3);

SimuParam.InvFov.inv_sc_pix_count = 501;          % 反演视场像素点个数
% SimuParam.InvFov.inv_sc_pix_count = 1001; 

%% 天线罩参数设置
% SimuParam.radome.thickness = 0.054; % 厚度
% SimuParam.radome.thickness = 0.15; % 厚度
SimuParam.radome.thickness = 0.06;   % 厚度
SimuParam.radome.epsilon = 3;       % 相对介电常数 
SimuParam.radome.tan_delta = 0.008; % 损耗角正切 调大增加衰减
SimuParam.radome.point_matrix = load('子系统2\参数预处理\参数\point_matrix.mat');  % 顶点矩阵
SimuParam.radome.patch_matrix = load('子系统2\参数预处理\参数\patch_matrix.mat');  % 面元矩阵
SimuParam.radome.T_radome_out = load('子系统2\参数预处理\参数\T_radome_out.mat');  % 温度分布
SimuParam.radome.T_radome_in = load('子系统2\参数预处理\参数\T_radome_in.mat'); 
SimuParam.radome.lay_num = 10; % 分层数

%% 输入场景设置
% SimuParam.Scene.fov_angle = SimuParam.Observe.vertical_scan * SimuParam.Observe.scan_step;  % 输入场景亮温图的视场范围
SimuParam.Scene.fov_angle = 10;

end
