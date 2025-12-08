% 该函数计算头罩的等效透过系数
% 输入：trans_h，phase_delay_h   水平极化电压透波率和插入相位延迟，以弧度形式传入
%       trans_v, phase_delay_v  垂直极化电压透波率和插入相位延迟
%       ant_polar   天线极化方向（天线接收电场方向）
%       n_radome    罩壁法向量
%       k_wave      波矢方向

% 输出：等效透过系数幅值trans_eq 和 相位phase_delay_eq

function [ trans_eq, phase_delay_eq ] = GetEquivalentTrans( trans_h, phase_delay_h, trans_v, phase_delay_v, ant_polar, n_radome, k_wave )
%% 计算入射平面法向量
n_in_plane = cross(n_radome, k_wave);

sin_alpha = dot(ant_polar,n_in_plane)/(sqrt(dot(ant_polar,ant_polar) * dot(n_in_plane,n_in_plane)));

cos_alpha = sqrt(1 - sin_alpha^2);

delta_phase_delay = phase_delay_h - phase_delay_v;

alpha_co = atan(trans_v*sin_alpha^2*sin(delta_phase_delay)/(trans_h*cos_alpha^2+trans_v*cos(delta_phase_delay)*sin_alpha^2));

trans_eq = sqrt(trans_h^2*cos_alpha^4 + trans_v^2*sin_alpha^4 + 2*trans_h*trans_v*cos(delta_phase_delay)*sin_alpha^2*cos_alpha^2);

phase_delay_eq = phase_delay_h - alpha_co;
