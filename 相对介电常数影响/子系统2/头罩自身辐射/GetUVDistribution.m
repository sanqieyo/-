function [ UV, UV_distribution, visib_redun_avg_matrix, redundant_baseline, trans_matrix] = GetUVDistribution( ant_pos, lambda )


%% 计算无冗余uv分布
ii=0;
uv = zeros(size(ant_pos,2),size(ant_pos,2)); % 全部的uv
for p = 1:size(ant_pos,2)
    for q = 1:size(ant_pos,2)
        ii=ii+1;
        delta_u = (ant_pos(1,p)-ant_pos(1,q))/lambda;
        delta_v = (ant_pos(2,p)-ant_pos(2,q))/lambda;
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


end

