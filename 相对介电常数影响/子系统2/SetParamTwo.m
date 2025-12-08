function SetParamTwo()
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
global SimuParam;

%% 参数预处理
[lambda, Gain, Omega] = ProcessSystemInput( SimuParam.SystemInput);
SimuParam.SystemInput.lambda = lambda;
SimuParam.SystemInput.Gain = Gain;
SimuParam.SystemInput.Omega = Omega;

[pattern_sc, ksai_pat_all, eta_pat_all, cosdtheta_pat_all, del_fov_pat, ksai_pat, eta_pat] = ProcessAntPattern( SimuParam.AntPattern);
SimuParam.AntPattern.pattern_sc = pattern_sc;
SimuParam.AntPattern.ksai_pat_all = ksai_pat_all;
SimuParam.AntPattern.eta_pat_all = eta_pat_all;
SimuParam.AntPattern.cosdtheta_pat_all = cosdtheta_pat_all;
SimuParam.AntPattern.del_fov_pat = del_fov_pat;
SimuParam.AntPattern.ksai_pat = ksai_pat;
SimuParam.AntPattern.eta_pat = eta_pat;

[del_u, del_v, del_uv, UV_distribution, UV, redundant_baseline, trans_matrix, array_num, Ant_array, ant_pos, visib_redun_avg_matrix] = ProcessArray( SimuParam.Array, SimuParam.SystemInput.lambda  );
SimuParam.Array.del_u = del_u;
SimuParam.Array.del_v = del_v;
SimuParam.Array.del_uv = del_uv;
SimuParam.Array.UV_distribution = UV_distribution;
SimuParam.Array.UV = UV;
SimuParam.Array.redundant_baseline = redundant_baseline;
SimuParam.Array.trans_matrix = trans_matrix;
SimuParam.Array.array_num = array_num;
SimuParam.Array.Ant_array = Ant_array;
SimuParam.Array.ant_pos = ant_pos;
SimuParam.Array.visib_redun_avg_matrix = visib_redun_avg_matrix;

[ksai_INV, eta_INV, ksai_inv, eta_inv, ksai_inv_rec, eta_inv_rec, in_blur] = ProcessInvFov(SimuParam.InvFov, SimuParam.Array.del_v);
SimuParam.InvFov.ksai_INV = ksai_INV;
SimuParam.InvFov.eta_INV = eta_INV;
SimuParam.InvFov.ksai_inv = ksai_inv;
SimuParam.InvFov.eta_inv = eta_inv;
SimuParam.InvFov.ksai_inv_rec = ksai_inv_rec;
SimuParam.InvFov.eta_inv_rec = eta_inv_rec;
SimuParam.InvFov.in_blur = in_blur;

[rank_point, normal_vector ] = ProcessRadome( SimuParam.radome );
SimuParam.radome.rank_point = rank_point;
SimuParam.radome.normal_vector = normal_vector;


end

