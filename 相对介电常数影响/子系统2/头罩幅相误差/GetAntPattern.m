function [ wait_cope_data,ksai_data,eta_data,ksai_pat_all, eta_pat_all, cosdtheta_pat_all, pattern_sc ,del_fov_pat ] = GetAntPattern( text, FOP_theta, pat_pix_count )
% 处理原始方向图数据，得到给定观测角和像素点数条件下的单元天线方向图

% 输入：     方向图数据 .txt文件
%            FOP_theta  天线观测角（正负60度时，这里就是60）
%            pat_pix_count 像素点数

% 输出：     ksai_data  即l = sin(theta)*cos(phi)
%            eta_data   即m = sin(theta)*sin(phi)
%            wait_cope_data 对应的功率增益

%            ksai_pat_all 给定观测角和像素点数条件下的l
%            eta_pat_all  给定观测角和像素点数条件下的m
%            cosdtheta_pat_all 给定观测角和像素点数条件下的倾斜因子
%            pattern_sc        给定观测角和像素点数条件下的归一化功率增益
%            del_fov_pat       给定观测角和像素点数条件下的像素点大小，即del_l * del_m

data = importdata(text);
%% 处理94new.txt原始方向图数据
wait_cope_data = [];
for i = 1 : size(data.data,1) / 721
    wait_cope_data = [wait_cope_data,data.data(1+(i-1)*721 : i*721,3)];
end   
phi_data = linspace(0,359.75,1440);
theta_data = linspace(0,180,721);
[theta_data,phi_data]=ndgrid(theta_data,phi_data);
ksai_data = sind(theta_data).*cosd(phi_data);
eta_data = sind(theta_data).*sind(phi_data);
%% 插值方向图
ksai_data = ksai_data(1 : 361, :);
eta_data = eta_data(1 : 361, :);
wait_cope_data = wait_cope_data(1 : 361, :);
wait_cope_data = 10.^(wait_cope_data/10);

%% 根据天线观测角度和像素点数确定匹配的方向图
Scene_Xi_pat = sind(FOP_theta);
ksai_pat=linspace(-Scene_Xi_pat,Scene_Xi_pat,pat_pix_count);
eta_pat =linspace(-Scene_Xi_pat,Scene_Xi_pat,pat_pix_count);
[ksai_pat_all,eta_pat_all]=ndgrid(ksai_pat,eta_pat);
del_fov_pat = (Scene_Xi_pat*2/length(ksai_pat)) * (Scene_Xi_pat*2/length(eta_pat));

cosdtheta_pat_all_temp = 1 - ksai_pat_all .* ksai_pat_all - eta_pat_all .* eta_pat_all;
cosdtheta_pat_all_temp(cosdtheta_pat_all_temp < 0) = nan; % 标记视场外点
cosdtheta_pat_all = real(sqrt(cosdtheta_pat_all_temp));
cosdtheta_pat_all(cosdtheta_pat_all < cosd(FOP_theta)) = nan;

pattern_sc = griddata(ksai_data,eta_data,wait_cope_data,ksai_pat_all,eta_pat_all);
pattern_sc = pattern_sc/max(max(pattern_sc)); % 归一化功率方向图
pattern_sc(isnan(pattern_sc)) = 0;



end

