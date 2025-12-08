%

% 本程序可用于完成复杂场景的亮温分布和天线温度的仿真计算
% 编写日期2016年5月14日
addpath(genpath(pwd))
clc;
clear all;
close all;
tic;
global SimuParam;


load('three_points.mat');


SetParam();                   % 设置仿真参数及参数预处理

% 子系统2运行函数
% 输入：    TB                   子系统1生成的亮温矩阵
% 输出：    Visib                头罩影响下的可见度

%           visibility_FF        头罩幅相误差影响的可见度
%           visibility_NF        头罩热辐射影响可见度
%           visibility_Ideal     没有头罩影响的理想可见度

[ ~, visibility_FF_1 , ~, visibility_Ideal_1] = SystemTwo( TB_1 );
[ ~, visibility_FF_2 , ~, visibility_Ideal_2] = SystemTwo( TB_2 );
[ ~, visibility_FF_3 , ~, visibility_Ideal_3] = SystemTwo( TB_3 );


% 对可见度反演
ksai_INV = SimuParam.InvFov.ksai_INV;
eta_INV = SimuParam.InvFov.eta_INV;
eta_inv= SimuParam.InvFov.eta_inv;
ksai_inv = SimuParam.InvFov.ksai_inv;
in_blur = SimuParam.InvFov.in_blur;
UV_distribution = SimuParam.Array.UV_distribution;
del_uv = SimuParam.Array.del_uv;
Omega = SimuParam.SystemInput.Omega;


TB_inv_FF_1 = TB_inv_process(ksai_INV,eta_INV,ksai_inv,eta_inv,in_blur,UV_distribution,del_uv,Omega,visibility_FF_1);
TB_inv_ideal_1 = TB_inv_process(ksai_INV,eta_INV,ksai_inv,eta_inv,in_blur,UV_distribution,del_uv,Omega,visibility_Ideal_1);

TB_inv_FF_2 = TB_inv_process(ksai_INV,eta_INV,ksai_inv,eta_inv,in_blur,UV_distribution,del_uv,Omega,visibility_FF_2);
TB_inv_ideal_2 = TB_inv_process(ksai_INV,eta_INV,ksai_inv,eta_inv,in_blur,UV_distribution,del_uv,Omega,visibility_Ideal_2);

TB_inv_FF_3 = TB_inv_process(ksai_INV,eta_INV,ksai_inv,eta_inv,in_blur,UV_distribution,del_uv,Omega,visibility_FF_3);
TB_inv_ideal_3 = TB_inv_process(ksai_INV,eta_INV,ksai_inv,eta_inv,in_blur,UV_distribution,del_uv,Omega,visibility_Ideal_3);




figure;
imagesc(ksai_INV, eta_INV, TB_inv_FF_1);
colormap(jet);
colorbar;
shading interp;
title('头罩影响的可见度--1');
figure;
imagesc(ksai_INV, eta_INV, TB_inv_ideal_1);
colormap(jet);
colorbar;
shading interp;
title('理想可见度--1');


figure;
imagesc(ksai_INV, eta_INV, TB_inv_FF_2);
colormap(jet);
% colorbar;
shading interp;
title('头罩影响的可见度--2');
figure;
imagesc(ksai_INV, eta_INV, TB_inv_ideal_2);
colormap(jet);
% colorbar;
shading interp;
title('理想可见度--2');


figure;
imagesc(ksai_INV, eta_INV, TB_inv_FF_3);
colormap(jet);
% colorbar;
shading interp;
title('头罩影响的可见度--3');
figure;
imagesc(ksai_INV, eta_INV, TB_inv_ideal_3);
colormap(jet);
% colorbar;
shading interp;
title('理想可见度--3');
    


%% 仿真参数
% for ii = 1 : 2
%
%
%     %% 亮温计算
%
%     [TB] = SystemOne(ii);           % 计算亮温分布
%     %     [ Visib,~, ~, ~] = SystemTwo( TB );  % 计算可见度输出
%
%     TBs{1, ii} = TB;
%
% end


toc