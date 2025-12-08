% 子系统2 包含头罩幅相误差的场景可见度计算代码 2022.04.15
function [visibility_FF] = GetFFVisibility( Param ,TB_original)

%% 输入场景参数处理
[ TB_scene, del_Fov_sc, ksai_sc, eta_sc, TB_origin, del_Fov_origin, ksai_origin, eta_origin, TB_error, ksai_error, eta_error] = ProcessTB( Param, TB_original);

%% 下一部分
Ant_array = Param.Array.Ant_array;
ant_polar = Param.AntPattern.polarization; % 接收天线极化方向
ant_pos = Param.Array.ant_pos;
UV_distribution = Param.Array.UV_distribution;
del_uv = Param.Array.del_uv;

f0 = Param.SystemInput.f0;
% 头罩参数
rank_point = Param.radome.rank_point;
normal_vector = Param.radome.normal_vector;
target_rank = Param.radome.patch_matrix;
target_point = Param.radome.point_matrix;
% thickness = Param.radome.thickness;
epsilon = Param.radome.epsilon;
tan_delta = Param.radome.tan_delta;
% d_flat = thickness;

%% 计算头罩等效透过系数
trans_eq = zeros(length(Ant_array),length(ksai_error));
phase_delay_eq = zeros(length(Ant_array),length(ksai_error));

% for i = 1:length(Ant_array)
%     inc_node = [ant_pos(i,:) 0]; % 入射方向矢量
%     
%     for j = 1:length(ksai_error)
%         
%         [ facet_list5, k_wave, ~ ] = GetIntersectionPoint( inc_node, ksai_error(j),eta_error(j), rank_point, normal_vector,target_rank.patch_matrix,target_point.point_matrix);
%         N_j = normal_vector(facet_list5,:);
%         theta_i = acos((dot(k_wave,N_j ))/(sqrt(dot(k_wave, k_wave) * dot(N_j, N_j))));
%         
%         if theta_i > pi/2
%             theta_i = pi - theta_i;
%         end
%         theta_i = theta_i/pi*180;
%         
%         patch_index = target_rank.patch_matrix(facet_list5,:);
%         if target_point.point_matrix(patch_index(1),3) > 0.8 || target_point.point_matrix(patch_index(2),3) > 0.8 || target_point.point_matrix(patch_index(3),3) > 0.8
%             thickness = 0.15;
%         else
%             thickness = Param.radome.thickness;
%         end
%         d_flat = thickness/cosd(theta_i);
% 
%         [ trans_h, phase_delay_h ] = GetHorizontalTrans( theta_i,d_flat,epsilon,tan_delta,f0);
%         [ trans_v, phase_delay_v ] = GetVerticalTrans( theta_i,d_flat,epsilon,tan_delta,f0);
%         [ trans_eq(i,j), phase_delay_eq(i,j) ] = GetEquivalentTrans( trans_h, phase_delay_h, trans_v, phase_delay_v, ant_polar, N_j, k_wave );
%     end
% end
% 
% save('trans_eq','trans_eq');
% save('phase_delay_eq','phase_delay_eq');

load('trans_eq');
load('phase_delay_eq');

mean_trans = (mean(mean(trans_eq)))^2;    % 平均功率传输系数

trans_eq_ori = zeros(length(Ant_array),length(ksai_origin));
phase_delay_eq_ori = zeros(length(Ant_array),length(ksai_origin));

Ksai_error = reshape(ksai_error,size(TB_error,1),size(TB_error,2));
Eta_error = reshape(eta_error,size(TB_error,1),size(TB_error,2));
Ksai_origin = reshape(ksai_origin,sqrt(length(TB_origin)),sqrt(length(TB_origin)));
Eta_origin = reshape(eta_origin,sqrt(length(TB_origin)),sqrt(length(TB_origin)));

parfor kk = 1:length(Ant_array)
    
    trans_eq_ori(kk,:) = reshape(griddata(Ksai_error,Eta_error,reshape(trans_eq(kk,:),size(TB_error,1),size(TB_error,2)),Ksai_origin, Eta_origin),1,[]); 

    phase_delay_eq_ori(kk,:) = reshape(griddata(Ksai_error,Eta_error,reshape(phase_delay_eq(kk,:),size(TB_error,1),size(TB_error,2)),Ksai_origin, Eta_origin),1,[]);

end


%% 获取无噪声可见度函数(加罩)
lambda = Param.SystemInput.lambda;
BandWidth = Param.SystemInput.BandWidth;
Omega = Param.SystemInput.Omega;
visib_redun_avg_matrix = Param.Array.visib_redun_avg_matrix;
UV = Param.Array.UV;
array_num = Param.Array.array_num;

Visibility_radome = zeros(length(Ant_array),length(Ant_array));
Visibility = zeros(length(Ant_array),length(Ant_array));
Visibility_large = zeros(length(Ant_array),length(Ant_array));
Visibility_large_radome = zeros(length(Ant_array),length(Ant_array));

% 可见度计算（小图＋大图）
for ii = 1 : length(Ant_array) 
    temp_pos = Ant_array(ii);
    temp_trans = trans_eq_ori(ii,:);
    temp_delay = phase_delay_eq_ori(ii,:);
    parfor jj = 1 : length(Ant_array)  
       u = real(temp_pos-Ant_array(jj))/lambda;
       v = imag(temp_pos-Ant_array(jj))/lambda;
       phi_delay = (temp_delay-phase_delay_eq_ori(jj,:))';
       fourier_vector_radome = del_Fov_origin.*exp(-1j*2*pi*(u*ksai_origin+v*eta_origin)-1j*phi_delay).*sinc(-BandWidth*(u*ksai_origin+v*eta_origin)/f0);
       Visibility_radome(ii,jj) = (temp_trans.*trans_eq_ori(jj,:).*TB_origin')*fourier_vector_radome/Omega; % 加罩可见度
       
%        fourier_vector = del_Fov_origin.*exp(-1j*2*pi*(u*ksai_origin+v*eta_origin)).*sinc(-BandWidth*(u*ksai_origin+v*eta_origin)/f0);
%        Visibility(i,j) = TB_origin'*fourier_vector/Omega;  % 未加罩可见度
       
       fourier_vector_large_radome = del_Fov_sc.*exp(-1j*2*pi*(u*ksai_sc+v*eta_sc)).*sinc(-BandWidth*(u*ksai_sc+v*eta_sc)/f0);
       Visibility_large_radome(ii,jj) = mean_trans.*TB_scene'*fourier_vector_large_radome/Omega;    % 有天线罩影响的背景可见度
       
%        fourier_vector_large = del_Fov_sc.*exp(-1j*2*pi*(u*ksai_sc+v*eta_sc)).*sinc(-BandWidth*(u*ksai_sc+v*eta_sc)/f0);
%        Visibility_large(i,j) = TB_scene'*fourier_vector_large/Omega;  % 没有天线罩影响的背景可见度
       
    end
end
% Visibility_radome = Visibility_radome-diag(diag(Visibility_radome)); % 加罩   减去零基线可见度
% Visibility = Visibility-diag(diag(Visibility)); % 未加罩 减去零基线可见度
% Visibility_large = Visibility_large-diag(diag(Visibility_large));

Visibility_radome = reshape(Visibility_radome, [], 1);
Visibility_large_radome = reshape(Visibility_large_radome,[],1);
% Visibility = reshape(Visibility, [], 1);
% Visibility_large = reshape(Visibility_large,[],1);


visibility_radome = visib_redun_avg_matrix*Visibility_radome; % 冗余平均
visibility_large_radome = visib_redun_avg_matrix*Visibility_large_radome;
% visibility = visib_redun_avg_matrix*Visibility; % 冗余平均
% visibility_large = visib_redun_avg_matrix*Visibility_large;

visibility_radome = visibility_radome + visibility_large_radome;  % 有天线罩影响的可见度
% visibility = visibility + visibility_large;         % 没有天线罩影响的可见度

visibility_FF = visibility_radome;

end

