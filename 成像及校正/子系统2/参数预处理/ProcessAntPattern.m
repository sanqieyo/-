function [ pattern_sc, ksai_pat_all, eta_pat_all, cosdtheta_pat_all, del_fov_pat, ksai_pat, eta_pat ] = ProcessAntPattern( Param )
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明

%% 导入方向图
data = importdata(Param.pat_file_path);  %% 94new.txt

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



FOP_theta = Param.FOP_theta;
cur_pat_pix_count = Param.pat_pix_count;
Scene_Xi_pat = sind(FOP_theta);
ksai_pat=linspace(-Scene_Xi_pat,Scene_Xi_pat,cur_pat_pix_count);
eta_pat =linspace(-Scene_Xi_pat,Scene_Xi_pat,cur_pat_pix_count);
[ksai_pat_all,eta_pat_all]=ndgrid(ksai_pat,eta_pat);
del_fov_pat = (Scene_Xi_pat*2/length(ksai_pat)) * (Scene_Xi_pat*2/length(eta_pat));

cosdtheta_pat_all_temp = 1 - ksai_pat_all .* ksai_pat_all - eta_pat_all .* eta_pat_all;
cosdtheta_pat_all_temp(cosdtheta_pat_all_temp < 0) = nan;%标记视场外点
cosdtheta_pat_all = real(sqrt(cosdtheta_pat_all_temp));
cosdtheta_pat_all(cosdtheta_pat_all < cosd(FOP_theta)) = nan;

pattern_sc = griddata(ksai_data,eta_data,wait_cope_data,ksai_pat_all,eta_pat_all);
pattern_sc = pattern_sc/max(max(pattern_sc));
pattern_sc(isnan(pattern_sc)) = 0;

end


