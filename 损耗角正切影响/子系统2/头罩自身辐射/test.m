clc
clear

%% 系统参数设置
c = 3.0e8;            % 光速       
f0 = 94e9;            % 中心频率
lambda = c/f0;        % 波长
% 系统噪声
Tau = 0.1;              % 积分时间
Ta = 265;             % 天线温度
Tr = 434;             % 接收机等效噪声温度
BandWidth = 2000.e6;   % 带宽，不考虑消条纹效应时可设为0
% gain = 25;            % 天线的功率增益/dB，方向性系数
gain = 23;            
Gain = 10^(gain/10);  % 天线功率倍数增益
Omega = 4*pi/Gain;    % 天线波束立体角
Tsys = Ta + Tr;       % 系统温度
del_vis = Tsys/sqrt(BandWidth*Tau);   % 测量可见度的标准差

%% 天线阵列，96通道
load('array.mat');    % 没有进行波长归一化
ant_pos = zeros(96,2);
ant_pos(:,1) = real(elem_pos_opt);
ant_pos(:,2) = imag(elem_pos_opt);
Ant_array = elem_pos_opt';
array_num = length(Ant_array);

del_u = 0.021/lambda;
del_v = 0.021*sqrt(3)/2/lambda;
del_uv = del_u*del_v;

%% 获取所有UV，无冗余基线UV_distribution，可见度冗余平均矩阵visib_redun_avg_matrix等
[ UV, UV_distribution, visib_redun_avg_matrix, redundant_baseline, trans_matrix] = GetUVDistribution( ant_pos', lambda );

%% 导入天线罩网格划分数据(几何参数)
load('radome_tri_partition_20mm.mat');
vertex_x = 10.^(-3).*radomepartition1(1:2:9679,3); % 所有顶点坐标
vertex_y = 10.^(-3).*radomepartition1(1:2:9679,4);
vertex_z = 10.^(-3).*radomepartition1(2:2:9680,3);
points_index = radomepartition1(9691:19366,4:6); % 三角面元个数，每一行代表三角面元的三个顶点索引

% figure
% scatter3(vertex_x,vertex_y,vertex_z,'*');
% axis equal

%% 所有三角面元的顶点坐标
triangular_patch_x = [vertex_x(points_index(:,1)) vertex_x(points_index(:,2)) vertex_x(points_index(:,3))]; % 9676 x 3
triangular_patch_y = [vertex_y(points_index(:,1)) vertex_y(points_index(:,2)) vertex_y(points_index(:,3))]; 
triangular_patch_z = [vertex_z(points_index(:,1)) vertex_z(points_index(:,2)) vertex_z(points_index(:,3))]; 

%% 所有三角面元中心点坐标
triangular_patch_center_x = sum(triangular_patch_x,2)./3;
triangular_patch_center_y = sum(triangular_patch_y,2)./3;
triangular_patch_center_z = sum(triangular_patch_z,2)./3;
% triangular_patch_center_pos = [triangular_patch_center_x triangular_patch_center_y triangular_patch_center_z];

patch_index = find(triangular_patch_center_z > 0.2);    % z > 0.2部分面元坐标
triangular_patch_center_x = triangular_patch_center_x(patch_index);
triangular_patch_center_y = triangular_patch_center_y(patch_index);
triangular_patch_center_z = triangular_patch_center_z(patch_index);
triangular_patch_center_pos = [triangular_patch_center_x triangular_patch_center_y triangular_patch_center_z];
triangular_patch_x = triangular_patch_x(patch_index,:);
triangular_patch_y = triangular_patch_y(patch_index,:);
triangular_patch_z = triangular_patch_z(patch_index,:);



patch_num = size(triangular_patch_x,1); %三角面元个数
z_direction_vector = [zeros(patch_num,1) zeros(patch_num,1) ones(patch_num,1)];
cos_theta_of_patch_z = triangular_patch_center_z./(sqrt(triangular_patch_center_x.^2 + triangular_patch_center_y.^2 ...
                                                        + triangular_patch_center_z.^2));                                                    
theta_deg = acosd(cos_theta_of_patch_z);  % 三角面元的天顶角  计算天线罩亮温时角度(不是)
% index = find(triangular_patch_center_pos(:,3) < 0.1);

% C = ones(size(triangular_patch_center_pos,1),1);  % 设置天线罩颜色，表示其亮温分布
% for i = 1:size(triangular_patch_center_pos,1)
%    if (triangular_patch_center_pos(i,3) >= 1)
%        C(i) = 800;
% %        elseif (triangular_patch_center_pos(i,3) >= 0.6)
% %            C(i) = 400;
% %            elseif (triangular_patch_center_pos(i,3) >= 0.2)
% %                C(i) = 300;
%                else
%                    C(i) = 0;
%    end
% end
% TB = C;
% figure  % 绘制天线罩网格划分三角面元
% for i = 1 : patch_num
%     patch(triangular_patch_x(i,:),triangular_patch_y(i,:),triangular_patch_z(i,:),C(i));
%     hold on
% end
% axis equal
% colormap(jet)

%% 三角面元法向量
edge_vector_1 = [triangular_patch_x(:,2)-triangular_patch_x(:,1) triangular_patch_y(:,2)-triangular_patch_y(:,1) triangular_patch_z(:,2)-triangular_patch_z(:,1)];% 三角面元边矢量
edge_vector_2 = [triangular_patch_x(:,3)-triangular_patch_x(:,1) triangular_patch_y(:,3)-triangular_patch_y(:,1) triangular_patch_z(:,3)-triangular_patch_z(:,1)];
edge_vector_3 = [triangular_patch_x(:,2)-triangular_patch_x(:,3) triangular_patch_y(:,2)-triangular_patch_y(:,3) triangular_patch_z(:,2)-triangular_patch_z(:,3)];
normal_vector = zeros(size(edge_vector_1,1), size(edge_vector_1,2));
for i = 1 : size(edge_vector_1,1)  % 叉积求三角面元法向量
    normal_vector(i,:) = cross(edge_vector_1(i,:), edge_vector_2(i,:));
end
norm = sqrt(normal_vector(:,1).^2 + normal_vector(:,2).^2 + normal_vector(:,3).^2);
normal_vector_norm = normal_vector./norm;  % 法向量归一化

%% 根据三角面元顶点计算面积
edge_length_1 = sqrt((edge_vector_1(:,1)).^2 + (edge_vector_1(:,2)).^2 + (edge_vector_1(:,3)).^2);
edge_length_2 = sqrt((edge_vector_2(:,1)).^2 + (edge_vector_2(:,2)).^2 + (edge_vector_2(:,3)).^2);
edge_length_3 = sqrt((edge_vector_3(:,1)).^2 + (edge_vector_3(:,2)).^2 + (edge_vector_3(:,3)).^2);
p = (edge_length_1+edge_length_2+edge_length_3)./2;
S_tri = sqrt(p.*(p-edge_length_1).*(p-edge_length_2).*(p-edge_length_3)); % 计算三角面元面积

%% 计算三角面元法向量 与 三角面元和天线连线 的夹角cos_theta
% cos_theta = zeros(patch_num, size(ant_pos,1));
% for i = 1 : patch_num
%     for j = 1 : size(ant_pos,1)
%         
%         vector_1 = [triangular_patch_center_pos(i,1)-ant_pos(j,1) triangular_patch_center_pos(i,2)-ant_pos(j,2) triangular_patch_center_pos(i,3)];
%         vector_2 = normal_vector_norm(i,:);
%         vector_1_length = sqrt(vector_1(1)^2 + vector_1(2)^2 + vector_1(3)^2);
%         vector_2_length = sqrt(vector_2(1)^2 + vector_2(2)^2 + vector_2(3)^2);
%         cos_theta(i,j) = (vector_1(1)*vector_2(1) + vector_1(2)*vector_2(2) + vector_1(3)*vector_2(3))/vector_1_length/vector_2_length;
%         
%     end
% end
% % theta_deg = acosd(cos_theta);


%% cos_theta应该这样计算：三角面元和天线连线 与 z轴 的夹角
%   天线看面元的方位角，天顶角
cos_theta = zeros(patch_num, size(ant_pos,1));  % 倾斜因子
cos_phi = zeros(patch_num, size(ant_pos,1));
for i = 1 : patch_num
    for j = 1 : size(ant_pos,1)
        
        vector_1 = [triangular_patch_center_pos(i,1)-ant_pos(j,1) triangular_patch_center_pos(i,2)-ant_pos(j,2) triangular_patch_center_pos(i,3)];
        vector_2 = [0 0 1];
        vector_1_length = sqrt(vector_1(1)^2 + vector_1(2)^2 + vector_1(3)^2);
        vector_2_length = sqrt(vector_2(1)^2 + vector_2(2)^2 + vector_2(3)^2);
        cos_theta(i,j) = (vector_1(1)*vector_2(1) + vector_1(2)*vector_2(2) + vector_1(3)*vector_2(3))/vector_1_length/vector_2_length;
        
        vector_1 = [1 - ant_pos(j,1) 0 0];
        vector_2 = [triangular_patch_center_pos(i,1)-ant_pos(j,1) triangular_patch_center_pos(i,2)-ant_pos(j,2) 0];
        vector_1_length = sqrt(vector_1(1)^2 + vector_1(2)^2 + vector_1(3)^2);
        vector_2_length = sqrt(vector_2(1)^2 + vector_2(2)^2 + vector_2(3)^2);
        cos_phi(i,j) = (vector_1(1)*vector_2(1) + vector_1(2)*vector_2(2) + vector_1(3)*vector_2(3))/vector_1_length/vector_2_length;
        
    end
end   

sin_theta = sqrt(1 - cos_theta.^2);
sin_phi = sqrt(1 - cos_phi.^2);
ksai_patch = sin_theta.*cos_phi;
eta_patch = sin_theta.*sin_phi;

%% 三角面元法向量 和 面元与原点连线 的夹角（新增，用于计算平板亮温的入射角）
cos_theta_in = zeros(patch_num, 1);
for i = 1 : patch_num
    
        vector_1 = normal_vector_norm(i,:);
        vector_2 = triangular_patch_center_pos(i,:);
        vector_1_length = sqrt(vector_1(1)^2 + vector_1(2)^2 + vector_1(3)^2);
        vector_2_length = sqrt(vector_2(1)^2 + vector_2(2)^2 + vector_2(3)^2);
        cos_theta_in(i) = (vector_1(1)*vector_2(1) + vector_1(2)*vector_2(2) + vector_1(3)*vector_2(3))/vector_1_length/vector_2_length;
        
end   % 入射角
theta_in = acos(cos_theta_in);

%% 三角面元法向量与z轴夹角(新增)
cos_theta_2 = zeros(patch_num, 1);
for i = 1 : patch_num        
        vector_1 = normal_vector(i,:);
        vector_2 = [0 0 1];
        vector_1_length = sqrt(vector_1(1)^2 + vector_1(2)^2 + vector_1(3)^2);
        vector_2_length = sqrt(vector_2(1)^2 + vector_2(2)^2 + vector_2(3)^2);
        cos_theta_2(i) = (vector_1(1)*vector_2(1) + vector_1(2)*vector_2(2) + vector_1(3)*vector_2(3))/vector_1_length/vector_2_length;
end
S_tri_new = S_tri.*abs(cos_theta_2);


%% 三角面元到各天线的距离   % 9676 x 24
distance = zeros(patch_num, size(ant_pos,1));
for i = 1 : patch_num
    for j = 1 : size(ant_pos,1)       
        distance(i,j) = sqrt((triangular_patch_center_pos(i,1)-ant_pos(j,1))^2 + (triangular_patch_center_pos(i,2)-ant_pos(j,2))^2 ...
                               + (triangular_patch_center_pos(i,3))^2);       
    end
end

%% 天线罩内外表面温度（温度分布）
T_radome_out = ones(size(triangular_patch_center_pos,1),1);  
T_radome_in = ones(size(triangular_patch_center_pos,1),1);  
for i = 1:size(triangular_patch_center_pos,1)
   if (triangular_patch_center_pos(i,3) >= 1.02)
       T_radome_out(i) = 1500*triangular_patch_center_pos(i,3);
       T_radome_in(i) = 500*triangular_patch_center_pos(i,3);
%        elseif (triangular_patch_center_pos(i,3) >= 0.6)
%            C(i) = 400;
%            elseif (triangular_patch_center_pos(i,3) >= 0.2)
%                C(i) = 300;
               else
                   T_radome_out(i) = 500*triangular_patch_center_pos(i,3);
                   T_radome_in(i) = 350*triangular_patch_center_pos(i,3);
   end
end
% TB = T_radome;

%% 计算天线罩表面反射率（材质参数）
c = 3.0e8;              
f0 = 94e9;          
lambda = c/f0;  
d_flat = 20*lambda; % 天线罩厚度
lay_num = 20;       % 介质分层数
delta_d = d_flat/lay_num;

epsilon = 3; 
delta = 0.003; 
epsilon_c = epsilon * (1 - 1j * delta);

alpha_2 = 2*pi/lambda*abs(imag(sqrt(epsilon_c)));
beta_2 = 2*pi/lambda*real(sqrt(epsilon_c));
k_alpha_2 = 2*alpha_2;  % 吸收系数

% 计算实透射角
k1 = 2*pi/lambda*sqrt(epsilon);
pp = 2*alpha_2*beta_2;
qq = beta_2^2 - alpha_2^2 - k1^2*(sin(theta_in)).^2;
theta_t = atan(sqrt(2)*k1*sin(theta_in)./(sqrt(sqrt(pp^2 + qq.^2) + qq)));

L2 = exp(k_alpha_2*delta_d./cos(theta_t));  % 损耗因子

% 计算功率反射系数，水平极化和垂直极化
Gamma_h = abs((cos(theta_in) - sqrt(epsilon_c - sin(theta_in).^2)) ./ (cos(theta_in) + sqrt(epsilon_c - sin(theta_in).^2))).^2;
Gamma_v = abs((epsilon_c * cos(theta_in) - sqrt(epsilon_c - sin(theta_in).^2)) ./ (epsilon_c * cos(theta_in) + sqrt(epsilon_c - sin(theta_in).^2))).^2;

%% 分别计算天线罩水平极化、垂直极化辐射亮温分布
T_layer = zeros(length(T_radome_out), lay_num);
Ts_layer = zeros(length(T_radome_out), lay_num+1);
N = lay_num+1;
TB_radome_h = zeros(length(T_radome_out) ,1);
TB_radome_v = zeros(length(T_radome_out) ,1);
for i = 1 : length(T_radome_out) % 面元个数循环
    T_layer(i,:) = GetProfileT( T_radome_out(i), T_radome_in(i), lay_num );
    Ts_layer(i,:) = [0 (1-1/L2(i))*T_layer(i,:)];
    
    % 计算每一层的辐射亮温
    Gamma = zeros(1,lay_num+1); % 反射率
    Gamma(1) = Gamma_h(i);
    Gamma(length(Gamma)) = Gamma_h(i);
    L = L2(i)*ones(1,lay_num+1);
    L(1) = 0;
    TB = zeros(1,N);
    for ii = 2 : N
        temp1 = 1;
        for j = 2 : ii
            temp1 = temp1 * (1-Gamma(j-1))/(1-Gamma(j-1)*Gamma(j)/L(j)/L(j));
        end

        temp2 = 1;
        for k = 2 : ii-1
            temp2 = temp2 / L(k);
        end

        TB(ii) = (1 + Gamma(ii)/L(ii)) * Ts_layer(i,ii) * temp1 * temp2;
    end
    TB_radome_h(i) = sum(TB);
end

for i = 1 : length(T_radome_out) % 面元个数循环
    T_layer(i,:) = GetProfileT( T_radome_out(i), T_radome_in(i), lay_num );
    Ts_layer(i,:) = [0 (1-1/L2(i))*T_layer(i,:)];
    
    % 计算每一层的辐射亮温
    Gamma = zeros(1,lay_num+1); % 反射率
    Gamma(1) = Gamma_v(i);
    Gamma(length(Gamma)) = Gamma_v(i);
    L = L2(i)*ones(1,lay_num+1);
    L(1) = 0;
    TB = zeros(1,N);
    for ii = 2 : N
        temp1 = 1;
        for j = 2 : ii
            temp1 = temp1 * (1-Gamma(j-1))/(1-Gamma(j-1)*Gamma(j)/L(j)/L(j));
        end

        temp2 = 1;
        for k = 2 : ii-1
            temp2 = temp2 / L(k);
        end

        TB(ii) = (1 + Gamma(ii)/L(ii)) * Ts_layer(i,ii) * temp1 * temp2;
    end
    TB_radome_v(i) = sum(TB);
end

%% 计算天线罩最终辐射亮温分布
% 计算极化角
ant_polar = [0 1 0]; % 天线极化方向ant_polar 面元法向量 normal_vector_norm
sin_alpha_p = zeros(size(normal_vector_norm,1), 1);
for i = 1 : size(normal_vector_norm,1)
    sin_alpha_p(i) = dot(ant_polar,normal_vector_norm(i,:))/(sqrt(dot(ant_polar,ant_polar) * dot(normal_vector_norm(i,:),normal_vector_norm(i,:))));
end
sin_alpha_p = abs(sin_alpha_p);
cos_alpha_p = sqrt(1 - sin_alpha_p.^2);
E_radome = sqrt(TB_radome_h).*cos_alpha_p + sqrt(TB_radome_v).*sin_alpha_p;
TB_radome = E_radome.^2;


%% 综合孔径近场可见度计算(未考虑天线方向图影响)
% Visibility = zeros(array_num,array_num);
% for i = 1 : array_num 
%     for j = 1 : array_num   
%         r1 = distance(:,i); 
%         r2 = distance(:,j);
% %         fourier_vector = S_tri.*sqrt(cos_theta(:,i).*cos_theta(:,j)).*exp(-1j*2*pi/lambda*(r1 - r2))./r1./r2;
%         fourier_vector = S_tri_new.*sqrt(cos_theta(:,i).*cos_theta(:,j)).*exp(-1j*2*pi/lambda*(r1 - r2))./r1./r2;
%         Visibility(i,j) = TB_radome'*fourier_vector/Omega;   
%     end
% end

%% 综合孔径近场可见度计算(考虑天线方向图影响)

%% 计算面元的theta、phi、l、m
% 面元位置坐标 triangular_patch_center_pos


%% 导入天线方向图
FOP_theta = 60;             % 天线方向图正负60度
sc_fov_angle = FOP_theta*2; % 场景观测角度（天线）
pat_pix_count = 500;        % 像素点500
text = '94new.txt';

[ wait_cope_data,ksai_data,eta_data,ksai_pat_all, eta_pat_all, cosdtheta_pat_all, pattern_sc, del_fov_pat ] = GetAntPattern( text, FOP_theta, pat_pix_count ); % 天线方向图

% figure % 绘制天线方向图
% surf(ksai_pat_all,eta_pat_all,pattern_sc);
% colorbar;
% colormap(jet);
% shading interp
% title('单元天线功率方向图');

ksai_patch_vec = reshape(ksai_patch, [], 1);
eta_patch_vec = reshape(eta_patch, [], 1);
pattern_patch_vec = griddata(ksai_pat_all,eta_pat_all,pattern_sc,ksai_patch_vec,eta_patch_vec); % 得到与TB相匹配的天线功率方向图
pattern_patch = reshape(pattern_patch_vec,size(ksai_patch,1),size(ksai_patch,2));
pattern_patch(isnan(pattern_patch)) = 0;


Visibility = zeros(array_num,array_num);
for i = 1 : array_num 
    for j = 1 : array_num   
        r1 = distance(:,i); 
        r2 = distance(:,j);
%         fourier_vector = S_tri.*sqrt(cos_theta(:,i).*cos_theta(:,j)).*exp(-1j*2*pi/lambda*(r1 - r2))./r1./r2;
        fourier_vector = S_tri_new.*sqrt(cos_theta(:,i).*cos_theta(:,j)).*exp(-1j*2*pi/lambda*(r1 - r2))./r1./r2;
        Visibility(i,j) = (sqrt(pattern_patch(:,i) .* pattern_patch(:,j)).*TB_radome)'*fourier_vector/Omega;    
    end
end

Visibility = Visibility - diag(diag(Visibility));  % 减了零基线怎么得到的结果亮温图数值更大？？
Visibility_v = reshape(Visibility,[],1);

visibility = visib_redun_avg_matrix*Visibility_v;
% visibility = Visibility_v;


%% 反演视场设置，获取反演参数
Scene_Xi_inv = sind(5)*2/sqrt(3); % 反演视场范围(矩形)
Scene_Eta_inv = sind(5)*2/sqrt(3);
% inv_sc_pix_count = 83; % 反演视场像素点个数
inv_sc_pix_count = 241;

% 综合孔径系统的分辨率，等于口径的两倍，输入场景的最小像素点应小于分辨率/主波束的1/3；
ksai_INV=linspace(-Scene_Xi_inv,Scene_Xi_inv,inv_sc_pix_count); % 场景的像素点大小必须小于系统主波束的1/3
eta_INV=linspace(-Scene_Eta_inv,Scene_Eta_inv,inv_sc_pix_count);% 此处修改像素数
[ksai_inv,eta_inv]=ndgrid(ksai_INV,eta_INV);
ksai_inv_rec = reshape(ksai_inv,[],1);
eta_inv_rec = reshape(eta_inv,[],1);

% del_u_new = del_u; % 这样画出来的六边形区域比较小
del_u_new = del_v; % 这样画出来的六边形区域比较大
x = linspace(0, 2*pi, 7);
X = exp(1j*x);
a = 1j*linspace(-1/del_u_new/2/sqrt(3),(1/del_u_new/2/sqrt(3)),20)+1/del_u_new/2;
fov2 = [a a*X(1) a*X(2) a*X(3) a*X(4) a*X(5) a*X(6)];
fov = fov2;
fov = fov*exp(1j*2*pi/12);
in_plot_figure = inpolygon(ksai_inv_rec,eta_inv_rec,real(fov),imag(fov));
in_blur = reshape(in_plot_figure,length(ksai_INV),length(eta_INV));


TB_inv = zeros(length(ksai_INV),length(eta_INV)); % 无冗余基线反演
for mm = 1 : size(ksai_inv,1)
    for nn = 1 : size(ksai_inv,2)
        Ksai = ksai_inv(mm,nn);
        Eta = eta_inv(mm,nn);
        ifourier_vector = exp(1j*2*pi*(real(UV_distribution)*Ksai+imag(UV_distribution)*Eta));
        TB_inv(mm,nn) = del_uv*Omega*visibility.'*ifourier_vector;
    end
end
TB_inv = real(TB_inv);

%% 让六边形视场外的TB值不显示
TB_inv = TB_inv.*in_blur;
[a b] = find(TB_inv==0);
for i = 1:length(a)
    TB_inv(a(i),b(i)) = nan;
end

figure
imagesc(ksai_INV, eta_INV, -TB_inv);
colormap(jet);
colorbar;
shading interp;

