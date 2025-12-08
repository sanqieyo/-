function [ del_u, del_v, del_uv, UV_distribution, UV, redundant_baseline, trans_matrix, array_num, Ant_array, ant_pos, visib_redun_avg_matrix] = ProcessArray( Param , lambda)
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明

array_mat = Param.array_mat;
load(array_mat);%实际名为elem_pos_opt

% lambda = SimuParam.SystemInput.lambda;
del_u_input = Param.del_u_input;
del_u = del_u_input/lambda;%最小UV
del_v = del_u * sqrt(3)/2;
del_uv = del_u*del_v;


Ant_array = elem_pos_opt.'; %导入的elem_pos_opt为一列  Ant_array为一行
ant_pos(:,1) = real(elem_pos_opt); %ant_pos两列
ant_pos(:,2) = imag(elem_pos_opt);

% Ant_array = elem_pos_opt;
% ant_pos(:,1) = real(elem_pos_opt.');
% ant_pos(:,2) = imag(elem_pos_opt.');

array_num = length(Ant_array);


%% 计算无冗余uv分布
Ant_pos = ant_pos';
ii=0;
uv = zeros(size(Ant_pos,2),size(Ant_pos,2)); % 全部的uv
for p = 1:size(Ant_pos,2)
    for q = 1:size(Ant_pos,2)
        ii=ii+1;
        delta_u = (Ant_pos(1,p)-Ant_pos(1,q))/lambda;
        delta_v = (Ant_pos(2,p)-Ant_pos(2,q))/lambda;
        uv(p,q) = delta_u + 1j*delta_v;
    end
end

UV = reshape(uv,[],1); 
UV_t = UV*10^6;
UV_T = round(UV_t);
[UV_Hex_t,N_index] = unique(UV_T, 'rows');
UV_distribution = UV_T(N_index,:);  % 获取无冗余uv：UV_distribution

%% 计算可见度冗余平均矩阵
redundant_baseline = zeros(length(UV_t),1);   % 基线平均矢量
trans_matrix = zeros(length(UV),length(UV_distribution));   % 零冗余基线向所有基线转换的矩阵，即 UV = trans_matrix*UV_distribution; 
for i = 1:length(N_index)
    [a b] = find(UV_T == UV_distribution(i));
    redundant_baseline(a) = 1/length(a);
    trans_matrix(a,i) = 1;
end
UV_distribution = UV_distribution/10^6;

redundant_avg_matrix = trans_matrix.*redundant_baseline;
visib_redun_avg_matrix = redundant_avg_matrix'; % 可见度冗余平均矩阵！！Visibility = visib_redun_avg_matrix' * R_matrix


