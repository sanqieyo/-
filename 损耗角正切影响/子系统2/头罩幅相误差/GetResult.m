function [visibility_radome, visibility, TB_inv_radome, TB_inv, amp_ratio, delta_theta] = GetResult(Param ,tan_delta)

% ##加入方向图影响，正演模型考虑14度视场##

%% 系统参数设置
c = 3.0e8;            % 光速       
f0 = 94e9;            % 中心频率
lambda = c/f0;        % 波长
% 系统噪声
Tau = 3;              % 积分时间
Ta = 100;             % 天线温度
Tr = 300;             % 接收机等效噪声温度
BandWidth = 400.e6;   % 带宽，不考虑消条纹效应时可设为0
gain = 23;            % 天线的功率增益/dB，方向性系数
Gain = 10^(gain/10);  % 天线功率倍数增益
Omega = 4*pi/Gain;    % 天线波束立体角
Tsys = Ta + Tr;       % 系统温度
del_vis = Tsys/sqrt(BandWidth*Tau);   % 测量可见度的标准差

%% 下一部分
Ant_array = Param.Array.Ant_array;
ant_polar = Param.AntPattern.polarization; % 接收天线极化方向
ant_pos = Param.Array.ant_pos;
UV_distribution = Param.Array.UV_distribution;
del_u = Param.Array.del_u;
del_v = Param.Array.del_v;
del_uv = Param.Array.del_uv;

f0 = Param.SystemInput.f0;
% 头罩参数
rank_point = Param.radome.rank_point;
normal_vector = Param.radome.normal_vector;
target_rank = Param.radome.patch_matrix;
target_point = Param.radome.point_matrix;
thickness = Param.radome.thickness;
epsilon = Param.radome.epsilon;
% tan_delta = Param.radome.tan_delta;
% d_flat = thickness;


%% 正演视场范围设置
% Scene_Xi_FOV = 1/del_v;  
% Scene_Eta_FOV = 1/del_v;
Scene_Xi_FOV = sind(5)*2/sqrt(3);
Scene_Eta_FOV = sind(5)*2/sqrt(3);
ksai_scene = linspace(-Scene_Xi_FOV,Scene_Xi_FOV,501);
eta_scene = linspace(-Scene_Eta_FOV,Scene_Eta_FOV,501);

[Ksai_scene,Eta_scene] = ndgrid(ksai_scene,eta_scene);
Ksai_scene_v = reshape(Ksai_scene,[],1);
Eta_scene_v = reshape(Eta_scene,[],1);    

del_ksai_scene = 2*Scene_Xi_FOV/length(ksai_scene);
del_eta_scene = 2*Scene_Eta_FOV/length(eta_scene);
del_Fov_sc = del_ksai_scene*del_eta_scene;  % 场景的像素的最小面积

%% 天线方向图
ksai_pat_all = Param.AntPattern.ksai_pat_all;
eta_pat_all = Param.AntPattern.eta_pat_all;
cosdtheta_pat_all = Param.AntPattern.cosdtheta_pat_all;
pattern_sc = Param.AntPattern.pattern_sc; % 只需要这个

%% 设置原始的输入场景
TB = 10*ones(length(ksai_scene),length(eta_scene));
TB(180,251) = 100000; 
Tb_v = reshape(TB,[],1);

%% 计算头罩等效透过系数
ksai_error = linspace(-Scene_Xi_FOV-0.01,Scene_Xi_FOV+0.01,51);
eta_error = linspace(-Scene_Xi_FOV-0.01,Scene_Xi_FOV+0.01,51);
[Ksai_error,Eta_error] = ndgrid(ksai_error,eta_error);
Ksai_error_v = reshape(Ksai_error,[],1);
 
trans_eq = zeros(length(Ant_array),length(Ksai_error_v));
phase_delay_eq = zeros(length(Ant_array),length(Ksai_error_v));

% thickness = d_radome;

for i = 1:length(Ant_array)
    inc_node = [ant_pos(i,:) 0]; % 入射方向矢量
    
    for j = 1:length(Ksai_error_v)
        
        [ facet_list5, k_wave, ~ ] = GetIntersectionPoint( inc_node, Ksai_error_v(j),Eta_scene_v(j), rank_point, normal_vector,target_rank.patch_matrix,target_point.point_matrix);
        N_j = normal_vector(facet_list5,:);
        theta_i = acos((dot(k_wave,N_j ))/(sqrt(dot(k_wave, k_wave) * dot(N_j, N_j))));
        
        if theta_i > pi/2
            theta_i = pi - theta_i;
        end
        theta_i = theta_i/pi*180;
        
        patch_index = target_rank.patch_matrix(facet_list5,:);
%         if target_point.point_matrix(patch_index(1),3) > 0.8 || target_point.point_matrix(patch_index(2),3) > 0.8 || target_point.point_matrix(patch_index(3),3) > 0.8
%             thickness = 0.15;
%         else
%             thickness = Param.radome.thickness;
%         end
        d_flat = thickness/cosd(theta_i);

        [ trans_h, phase_delay_h ] = GetHorizontalTrans( theta_i,d_flat,epsilon,tan_delta,f0);
        [ trans_v, phase_delay_v ] = GetVerticalTrans( theta_i,d_flat,epsilon,tan_delta,f0);
        [ trans_eq(i,j), phase_delay_eq(i,j) ] = GetEquivalentTrans( trans_h, phase_delay_h, trans_v, phase_delay_v, ant_polar, N_j, k_wave );
    end
end

save('trans_eq','trans_eq');
save('phase_delay_eq','phase_delay_eq');
% load('D:\我的学业\XH\02_DF26头罩影响仿真\trans_eq');
% load('D:\我的学业\XH\02_DF26头罩影响仿真\phase_delay_eq');

trans_eq_ori = zeros(length(Ant_array),length(Ksai_scene_v));
phase_delay_eq_ori = zeros(length(Ant_array),length(Ksai_scene_v));

parfor kk = 1:length(Ant_array)
    
    trans_eq_ori(kk,:) = reshape(griddata(Ksai_error,Eta_error,reshape(trans_eq(kk,:),size(Ksai_error,1),size(Ksai_error,2)),Ksai_scene, Eta_scene),1,[]); 

    phase_delay_eq_ori(kk,:) = reshape(griddata(Ksai_error,Eta_error,reshape(phase_delay_eq(kk,:),size(Ksai_error,1),size(Ksai_error,2)),Ksai_scene, Eta_scene),1,[]);

end

%% 分别计算可见度
lambda = Param.SystemInput.lambda;
BandWidth = Param.SystemInput.BandWidth;
Omega = Param.SystemInput.Omega;
visib_redun_avg_matrix = Param.Array.visib_redun_avg_matrix;
UV = Param.Array.UV;
array_num = Param.Array.array_num;

Visibility_radome = zeros(length(Ant_array),length(Ant_array));
Visibility = zeros(length(Ant_array),length(Ant_array));
for ii = 1 : length(Ant_array) 
    temp_pos = Ant_array(ii);
    temp_trans = trans_eq_ori(ii,:);
    temp_delay = phase_delay_eq_ori(ii,:);
    parfor jj = 1 : length(Ant_array)  
       u = real(temp_pos-Ant_array(jj))/lambda;
       v = imag(temp_pos-Ant_array(jj))/lambda;
       phi_delay = (temp_delay-phase_delay_eq_ori(jj,:))';
       fourier_vector_radome = del_Fov_sc.*exp(-1j*2*pi*(u*Ksai_scene_v+v*Eta_scene_v)+1j*phi_delay).*sinc(-BandWidth*(u*Ksai_scene_v+v*Eta_scene_v)/f0);
       Visibility_radome(ii,jj) = (temp_trans.*trans_eq_ori(jj,:).*Tb_v')*fourier_vector_radome/Omega; % 加罩可见度

       fourier_vector = del_Fov_sc.*exp(-1j*2*pi*(u*Ksai_scene_v+v*Eta_scene_v)).*sinc(-BandWidth*(u*Ksai_scene_v+v*Eta_scene_v)/f0);
       Visibility(ii,jj) = Tb_v'*fourier_vector/Omega;  % 未加罩可见度
              
    end
end

Visibility_radome = reshape(Visibility_radome, [], 1);
Visibility = reshape(Visibility,[],1);

visibility_radome = visib_redun_avg_matrix*Visibility_radome; % 冗余平均
visibility = visib_redun_avg_matrix*Visibility;

%% 反演视场范围

ksai_INV = Param.InvFov.ksai_INV;
eta_INV = Param.InvFov.eta_INV;
eta_inv= Param.InvFov.eta_inv;
ksai_inv = Param.InvFov.ksai_inv;
in_blur = Param.InvFov.in_blur;

TB_inv = zeros(length(ksai_INV),length(eta_INV)); % 无冗余基线反演
TB_inv_radome = zeros(length(ksai_INV),length(eta_INV));

for mm = 1 : size(ksai_inv,1)
    for nn = 1 : size(ksai_inv,2)
        Ksai = ksai_inv(mm,nn);
        Eta = eta_inv(mm,nn);
        ifourier_vector = exp(1j*2*pi*(real(UV_distribution)*Ksai+imag(UV_distribution)*Eta));
        TB_inv(mm,nn) = del_uv*Omega*visibility.'*ifourier_vector;
        TB_inv_radome(mm,nn) = del_uv*Omega*visibility_radome.'*ifourier_vector;
    end
end
TB_inv = real(TB_inv);
TB_inv_radome = real(TB_inv_radome);

% 让六边形视场外的TB值不显示
TB_inv = TB_inv.*in_blur;
TB_inv_radome = TB_inv_radome.*in_blur;
[a, b] = find(TB_inv==0);
for i = 1:length(a)
    TB_inv(a(i),b(i)) = nan;
end
[a, b] = find(TB_inv_radome==0);
for i = 1:length(a)
    TB_inv_radome(a(i),b(i)) = nan;
end

% figure
% imagesc(ksai_INV,eta_INV, TB_inv_radome);
% colormap(jet);
% colorbar;
% shading interp;
% 
% figure
% imagesc(ksai_INV,eta_INV, TB_inv);
% colormap(jet);
% colorbar;
% shading interp;

%% 计算一下衰减和位置偏移
max_TB = max(max(TB_inv));
max_TB_radome = max(max(TB_inv_radome));
amp_ratio = max_TB_radome/max_TB; % 衰减
[index_a,index_b] = find(TB_inv == max_TB);
eta_of_max_TB = eta_inv(index_a,index_b);
ksai_of_max_TB = ksai_inv(index_a,index_b);
theta_1 = asind(sqrt(eta_of_max_TB^2 + ksai_of_max_TB^2));
[index_a,index_b] = find(TB_inv_radome == max_TB_radome);
eta_of_max_TB_radome = eta_inv(index_a,index_b);
ksai_of_max_TB_radome = ksai_inv(index_a,index_b);
theta_2 = asind(sqrt(eta_of_max_TB_radome^2 + ksai_of_max_TB_radome^2));
delta_theta = abs(theta_2 - theta_1); %位置偏移



end

