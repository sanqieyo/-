% 子系统2 头罩自身辐射近场可见度计算代码 2022.03.07
function [ NFVisibility] = GetNFVisibility( Param )


points_index = Param.radome.patch_matrix;
point_index = Param.radome.point_matrix;
vertex_x = point_index.point_matrix(:,1);
vertex_y = point_index.point_matrix(:,2);
vertex_z = point_index.point_matrix(:,3);

%% 所有三角面元的顶点坐标
triangular_patch_x = [vertex_x(points_index.patch_matrix(:,1)) vertex_x(points_index.patch_matrix(:,2)) vertex_x(points_index.patch_matrix(:,3))]; % 9676 x 3
triangular_patch_y = [vertex_y(points_index.patch_matrix(:,1)) vertex_y(points_index.patch_matrix(:,2)) vertex_y(points_index.patch_matrix(:,3))]; 
triangular_patch_z = [vertex_z(points_index.patch_matrix(:,1)) vertex_z(points_index.patch_matrix(:,2)) vertex_z(points_index.patch_matrix(:,3))]; 

%% 所有三角面元中心点坐标
triangular_patch_center_x = sum(triangular_patch_x,2)./3;
triangular_patch_center_y = sum(triangular_patch_y,2)./3;
triangular_patch_center_z = sum(triangular_patch_z,2)./3;
% triangular_patch_center_pos = [triangular_patch_center_x triangular_patch_center_y triangular_patch_center_z];

%% 天线罩内外表面温度（温度分布）
T_radome_out = Param.radome.T_radome_out;  
T_radome_in = Param.radome.T_radome_in;  
normal_vector = Param.radome.normal_vector;

patch_index = find(triangular_patch_center_z > 0.001);    % z > 0.001部分面元坐标
triangular_patch_center_x = triangular_patch_center_x(patch_index);
triangular_patch_center_y = triangular_patch_center_y(patch_index);
triangular_patch_center_z = triangular_patch_center_z(patch_index);
triangular_patch_center_pos = [triangular_patch_center_x triangular_patch_center_y triangular_patch_center_z];
triangular_patch_x = triangular_patch_x(patch_index,:);
triangular_patch_y = triangular_patch_y(patch_index,:);
triangular_patch_z = triangular_patch_z(patch_index,:);
T_radome_out = T_radome_out.T_out(patch_index);
T_radome_in = T_radome_in.T_in(patch_index);
normal_vector = normal_vector(patch_index,:);


patch_num = size(triangular_patch_x,1); %三角面元个数
z_direction_vector = [zeros(patch_num,1) zeros(patch_num,1) ones(patch_num,1)];
cos_theta_of_patch_z = triangular_patch_center_z./(sqrt(triangular_patch_center_x.^2 + triangular_patch_center_y.^2 ...
                                                        + triangular_patch_center_z.^2));                                                    
theta_deg = acosd(cos_theta_of_patch_z);  % 三角面元的天顶角

%% 三角面元法向量
edge_vector_1 = [triangular_patch_x(:,2)-triangular_patch_x(:,1) triangular_patch_y(:,2)-triangular_patch_y(:,1) triangular_patch_z(:,2)-triangular_patch_z(:,1)];% 三角面元边矢量
edge_vector_2 = [triangular_patch_x(:,3)-triangular_patch_x(:,1) triangular_patch_y(:,3)-triangular_patch_y(:,1) triangular_patch_z(:,3)-triangular_patch_z(:,1)];
edge_vector_3 = [triangular_patch_x(:,2)-triangular_patch_x(:,3) triangular_patch_y(:,2)-triangular_patch_y(:,3) triangular_patch_z(:,2)-triangular_patch_z(:,3)];
% normal_vector = zeros(size(edge_vector_1,1), size(edge_vector_1,2));
% for i = 1 : size(edge_vector_1,1)  % 叉积求三角面元法向量
%     normal_vector(i,:) = cross(edge_vector_1(i,:), edge_vector_2(i,:));
% end
% norm = sqrt(normal_vector(:,1).^2 + normal_vector(:,2).^2 + normal_vector(:,3).^2);
% normal_vector_norm = normal_vector./norm;  % 法向量归一化

%% 根据三角面元顶点计算面积
edge_length_1 = sqrt((edge_vector_1(:,1)).^2 + (edge_vector_1(:,2)).^2 + (edge_vector_1(:,3)).^2);
edge_length_2 = sqrt((edge_vector_2(:,1)).^2 + (edge_vector_2(:,2)).^2 + (edge_vector_2(:,3)).^2);
edge_length_3 = sqrt((edge_vector_3(:,1)).^2 + (edge_vector_3(:,2)).^2 + (edge_vector_3(:,3)).^2);
p = (edge_length_1+edge_length_2+edge_length_3)./2;
S_tri = sqrt(p.*(p-edge_length_1).*(p-edge_length_2).*(p-edge_length_3)); % 计算三角面元面积

%% cos_theta应该这样计算：三角面元和天线连线 与 z轴 的夹角
%   天线看面元的方位角，天顶角
ant_pos = Param.Array.ant_pos;
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
    
        vector_1 = normal_vector(i,:);
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


%% 计算天线罩表面反射率（材质参数）
[ Gamma_h, Gamma_v, L2, ~ ] = GetGamma( Param, theta_in);

%% 分别计算天线罩水平极化、垂直极化辐射亮温分布
lay_num = Param.radome.lay_num;
T_layer = zeros(length(T_radome_out), lay_num);
Ts_layer = zeros(length(T_radome_out), lay_num+1); % 每一层的辐射能量
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
% ant_polar = [0 1 0];   % 天线极化方向ant_polar 面元法向量 normal_vector_norm
ant_polar = Param.AntPattern.polarization;
sin_alpha_p = zeros(size(normal_vector,1), 1);
for i = 1 : size(normal_vector,1)
    sin_alpha_p(i) = dot(ant_polar,normal_vector(i,:))/(sqrt(dot(ant_polar,ant_polar) * dot(normal_vector(i,:),normal_vector(i,:))));
end
sin_alpha_p = abs(sin_alpha_p);
cos_alpha_p = sqrt(1 - sin_alpha_p.^2);
E_radome = sqrt(TB_radome_h).*cos_alpha_p + sqrt(TB_radome_v).*sin_alpha_p;
TB_radome = E_radome.^2;



%% 综合孔径近场可见度计算(考虑天线方向图影响)

%% 导入天线方向图
ksai_pat_all = Param.AntPattern.ksai_pat_all;
eta_pat_all = Param.AntPattern.eta_pat_all;
pattern_sc = Param.AntPattern.pattern_sc;

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


lambda = Param.SystemInput.lambda;
array_num = Param.Array.array_num;
BandWidth = Param.SystemInput.BandWidth;
f0 = Param.SystemInput.f0;
Omega = Param.SystemInput.Omega;
visib_redun_avg_matrix = Param.Array.visib_redun_avg_matrix;
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

% Visibility = Visibility - diag(diag(Visibility));      
Visibility_v = reshape(Visibility,[],1);
visibility = visib_redun_avg_matrix*Visibility_v;
% visibility = Visibility_v;

NFVisibility = visibility;


end