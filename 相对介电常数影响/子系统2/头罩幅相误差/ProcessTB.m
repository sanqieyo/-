function [ TB_scene, del_Fov_sc, ksai_sc, eta_sc, TB_origin, del_Fov_origin, ksai_origin, eta_origin, TB_error,ksai_error, eta_error] = ProcessTB( Param, TB_original)
% 处理输入亮温数据，得到以阵列为中心的亮温矩阵以及均匀大背景矩阵
% 给定亮温矩阵TB和天线观测角度，得到相应的l，m以及倾斜因子

% 输入：     TB_original  子系统1输出的亮温矩阵

% 输出：     TB_scene     大背景亮温矩阵
%            del_Fov_sc   像素点大小，即del_l * del_m
%            ksai_sc      即l = sin(theta)*cos(phi)
%            eta_sc       即m = sin(theta)*sin(phi)

%            TB_origin         原始亮温图-均匀背景
%            del_Fov_origin    像素点大小，即del_l * del_m
%            ksai_origin       即l = sin(theta)*cos(phi)
%            eta_origin        即m = sin(theta)*sin(phi)

FOP_theta = Param.AntPattern.FOP_theta;
tar_fov_angle = FOP_theta*2; % 天线方向图角度120度
tar_fov_original = Param.Scene.fov_angle; % 正负5度原始视场

[ TB_large,TB_original,TB_error] = ProcessTBStart( Param, TB_original);  % 大背景图 和 原始图-均匀背景
% 分别计算两张图的正演参数
[tar_ksai, tar_eta, tar_offset, ~, ~] = GetOriginalLM( TB_large,tar_fov_angle );
[original_ksai, original_eta, original_offset, ~, ~] = GetOriginalLM( TB_original,tar_fov_original );% TB_tar主要是给出输出参数的size
[error_ksai, error_eta, error_offset, ~, ~] = GetOriginalLM( TB_error,tar_fov_original+1 ); % 用于计算天线罩误差

% 将L，M和倾斜因子转换为阵列为中心坐标系中
y_rotate_angle = 90; % 实际是逆时针旋转90度
y_rotate_matrix = [cosd(y_rotate_angle),sind(y_rotate_angle),0;
                 -sind(y_rotate_angle),cosd(y_rotate_angle), 0;
                 0,0,1];      % 绕y轴顺时针旋转矩阵，对应坐标为(z,x,y)             
z_rotate_angle = 90; % 实际是逆时针旋转90度
z_rotate_matrix = [1,0,0;
                    0,cosd(z_rotate_angle), sind(z_rotate_angle);
                    0,-sind(z_rotate_angle), cosd(z_rotate_angle)];% 绕y轴逆时针旋转矩阵，对应坐标为(z,x,y)
rotate_matrix = z_rotate_matrix * y_rotate_matrix;

tar_coefficeint = rotate_matrix * [reshape(tar_offset,1,[]);reshape(tar_ksai,1,[]);reshape(tar_eta,1,[])]; % 对应的顺序为sqrt(1-l^2-m^2)、l、m
ksai_tar_array_center = reshape(tar_coefficeint(2,:),size(tar_ksai)); % 原始TB对应的l，m和sqrt(1-l^2-m^2)坐标系旋转结果
eta_tar_array_center = reshape(tar_coefficeint(3,:),size(tar_eta));          
% cosdtheta_tar_array_center = reshape(tar_coefficeint(1,:),size(tar_eta));

tar_coefficeint_2 = rotate_matrix * [reshape(original_offset,1,[]);reshape(original_ksai,1,[]);reshape(original_eta,1,[])]; % 对应的顺序为sqrt(1-l^2-m^2)、l、m
ksai_tar_array_center_original = reshape(tar_coefficeint_2(2,:),size(original_ksai)); % 原始TB对应的l，m和sqrt(1-l^2-m^2)坐标系旋转结果
eta_tar_array_center_original = reshape(tar_coefficeint_2(3,:),size(original_eta)); 

tar_coefficeint_3 = rotate_matrix * [reshape(error_offset,1,[]);reshape(error_ksai,1,[]);reshape(error_eta,1,[])]; % 对应的顺序为sqrt(1-l^2-m^2)、l、m
ksai_tar_array_center_error = reshape(tar_coefficeint_3(2,:),size(error_ksai)); % 原始TB对应的l，m和sqrt(1-l^2-m^2)坐标系旋转结果
eta_tar_array_center_error = reshape(tar_coefficeint_3(3,:),size(error_eta)); 

%% 以TB_tar对应的l、m、del_fov确定相应的天线方向图（让天线方向图适应包含目标的小背景亮温图）
ksai_pat_all = Param.AntPattern.ksai_pat_all;
eta_pat_all = Param.AntPattern.eta_pat_all;
cosdtheta_pat_all = Param.AntPattern.cosdtheta_pat_all;
pattern_sc = Param.AntPattern.pattern_sc;

pattern_tar = griddata(ksai_pat_all,eta_pat_all,pattern_sc,ksai_tar_array_center,eta_tar_array_center); % 得到与TB相匹配的天线方向图
cosdtheta_tar = griddata(ksai_pat_all,eta_pat_all,cosdtheta_pat_all,ksai_tar_array_center,eta_tar_array_center); % 得到与TB相匹配的天线方向图
% TB_target = TB_original.*pattern_tar./ cosdtheta_tar;  % 经过方向图衰减的TB
TB_target = TB_large.*pattern_tar./ cosdtheta_tar;  % 经过方向图衰减的TB
TB_target(isnan(TB_target)) = 0;
TB_target_v = reshape(TB_target,[],1);

pattern_tar = griddata(ksai_pat_all,eta_pat_all,pattern_sc,ksai_tar_array_center_original,eta_tar_array_center_original); % 得到与TB相匹配的天线方向图
cosdtheta_tar = griddata(ksai_pat_all,eta_pat_all,cosdtheta_pat_all,ksai_tar_array_center_original,eta_tar_array_center_original); % 得到与TB相匹配的天线方向图
% TB_target = TB_original.*pattern_tar./ cosdtheta_tar;  % 经过方向图衰减的TB
TB_target_original = TB_original.*pattern_tar./ cosdtheta_tar;  % 经过方向图衰减的TB
TB_target_original(isnan(TB_target_original)) = 0;
TB_target_original_v = reshape(TB_target_original,[],1);

kernel_ksai = [0,0,0;
            -1,1,0;
            0,0,0];
kernel_eta = [0,-1,0;
            0,1,0;
            0,0,0];
co_x1 = conv2(ksai_tar_array_center,kernel_ksai,'same');
co_x1(:,size(co_x1,2)) = co_x1(:,size(co_x1,2)-1);
co_y1 = conv2(eta_tar_array_center,kernel_eta,'same');
co_y1(size(co_y1,1),:) = co_y1(size(co_y1,1)-1,:);
del_Fov_tar = abs(co_x1.*co_y1); 
ksai_tar_array_center_v = reshape(ksai_tar_array_center,[],1);
eta_tar_array_center_v = reshape(eta_tar_array_center,[],1);
del_Fov_tar_v = reshape(del_Fov_tar,[],1);    % 不均匀FOV

TB_scene = TB_target_v;
del_Fov_sc = del_Fov_tar_v;
ksai_sc = ksai_tar_array_center_v;
eta_sc = eta_tar_array_center_v;

co_x1 = conv2(ksai_tar_array_center_original,kernel_ksai,'same');
co_x1(:,size(co_x1,2)) = co_x1(:,size(co_x1,2)-1);
co_y1 = conv2(eta_tar_array_center_original,kernel_eta,'same');
co_y1(size(co_y1,1),:) = co_y1(size(co_y1,1)-1,:);
del_Fov_tar_original = abs(co_x1.*co_y1); 
ksai_tar_array_center_v = reshape(ksai_tar_array_center_original,[],1);
eta_tar_array_center_v = reshape(eta_tar_array_center_original,[],1);
del_Fov_tar_v = reshape(del_Fov_tar_original,[],1);    % 不均匀FOV

TB_origin = TB_target_original_v;
del_Fov_origin = del_Fov_tar_v;
ksai_origin = ksai_tar_array_center_v;
eta_origin = eta_tar_array_center_v;
% 上面为止部分是针对一张图的正确处理

ksai_tar_array_center_error_v = reshape(ksai_tar_array_center_error,[],1);
eta_tar_array_center_error_v = reshape(eta_tar_array_center_error,[],1);

ksai_error = ksai_tar_array_center_error_v;
eta_error = eta_tar_array_center_error_v;



%% 把正演参数都插值成均匀的
% tar_fov = Param.Scene.fov_angle;
% ksai_tar_new = linspace(-sind(tar_fov/2), sind(tar_fov/2), 101);
% eta_tar_new = linspace(-sind(tar_fov/2), sind(tar_fov/2), 101);
% [Ksai_tar_new, Eta_tar_new] = ndgrid(ksai_tar_new, eta_tar_new);
% ksai_tar_v = reshape(Ksai_tar_new, [] ,1);
% eta_tar_v = reshape(Eta_tar_new, [] ,1);
% del_ksai = 2*sind(tar_fov/2)/length(ksai_tar_new);
% del_eta = 2*sind(tar_fov/2)/length(eta_tar_new);
% del_Fov_tar = del_ksai*del_eta;      % 场景的像素的最小面积
% 
% TB_target_new = griddata(ksai_tar_array_center,eta_tar_array_center,TB_target,Ksai_tar_new, Eta_tar_new); 
% TB_target_new = reshape(TB_target_new,[],1);
% a = reshape(TB_target_new,length(ksai_tar_new),[]);


end

