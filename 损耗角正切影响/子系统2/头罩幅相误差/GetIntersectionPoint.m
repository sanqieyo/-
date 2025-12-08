function [ facet_list5, k_wave,source_point ] = GetIntersectionPoint( inc_node, l0, m0,rank_point, normal_vector,target_rank,target_point)
% 求三角面元与射线交点


z1 = 2;
r = z1/sqrt(1-l0^2-m0^2);
x1 = r*l0;
y1 = r*m0;

% x1 = L0*l0;
% y1 = L0*m0;
% z1 = L0*sqrt(1-l0^2-m0^2);
% 确定波矢
R_P0_P1 = sqrt((inc_node(1)-x1)^2+(inc_node(2)-y1)^2+z1^2);
kx = (x1-inc_node(1))/R_P0_P1;
ky = (y1-inc_node(2))/R_P0_P1;
kz = z1/R_P0_P1;
k_wave = [kx ky kz];


% R_source = 1;
source_point = 10*[l0 m0 sqrt(1-l0^2-m0^2)];
% 确定波矢
% k_wave = source_point - inc_node;

% inc_node = [ant_pos 0];
[facet_list5] = FindFacet2(k_wave,inc_node,rank_point,normal_vector,target_rank,target_point);

% if dot(k_wave,N_j) == 0  
%     P_j0 = [0 0 0];
% else
%     alpha = (dot(N_j,P_j1) - dot(N_j,P_0))/dot(k_wave,N_j);
%     P_j0 = P_0 + alpha * k_wave;
%     % 验证交点的合法性
%     % 计算各边的长度
%     P_j1_P_j2 = normest(P_j1 - P_j2);
%     P_j1_P_j3 = normest(P_j1 - P_j3);
%     P_j2_P_j3 = normest(P_j2 - P_j3);
%     P_j0_P_j1 = normest(P_j0 - P_j1);
%     P_j0_P_j2 = normest(P_j0 - P_j2);
%     P_j0_P_j3 = normest(P_j0 - P_j3);
%     
%     p1 = (P_j0_P_j1 + P_j0_P_j3 + P_j1_P_j3)/2;
%     p2 = (P_j0_P_j2 + P_j0_P_j3 + P_j2_P_j3)/2;
%     p3 = (P_j0_P_j1 + P_j0_P_j2 + P_j1_P_j2)/2;
%     p4 = (P_j1_P_j2 + P_j1_P_j3 + P_j2_P_j3)/2;
%     
%     S_1 = sqrt(p1*(p1-P_j0_P_j1)*(p1-P_j0_P_j3)*(p1-P_j1_P_j3));
%     S_2 = sqrt(p2*(p2-P_j0_P_j2)*(p2-P_j0_P_j3)*(p2-P_j2_P_j3));
%     S_3 = sqrt(p3*(p3-P_j0_P_j1)*(p3-P_j0_P_j2)*(p3-P_j1_P_j2));
%     S_4 = sqrt(p4*(p4-P_j1_P_j2)*(p4-P_j1_P_j3)*(p4-P_j2_P_j3));
%     
%     if S_1+S_2+S_3 ~= S_4
%         P_j0 = [0 0 0];
%     end
% end


% k_wave 与 N_j 夹角为 theta_i
        
end
