 clc
clear
%

% 本程序可用于完成复杂场景的亮温分布和天线温度的仿真计算
% 编写日期2016年5月14日
addpath(genpath(pwd))
clc;
clear all;
close all;
tic;
global SimuParam;

% 亮温场景输入

SetParam();                   % 设置仿真参数及参数预处理
addpath(genpath(pwd))
SetParamTwo();

d_radome = 0.06;
[visibility_radome, visibility, TB_inv_radome, TB_inv, amp_ratio, delta_theta] = GetResult(SimuParam,d_radome);

% 厚度影响
% d = 0:0.005:0.1;
% 
% amplitude_ratio = zeros(1,length(d));
% deviation_theta = zeros(1,length(d));
% parfor kk = 1 : length(d)
%     [~, ~, ~, ~,  amp_ratio, delta_theta ] = GetResult(SimuParam,d(kk))
%     amplitude_ratio(kk) = amp_ratio;
%     deviation_theta(kk) = delta_theta;
% end
% 
% save('amplitude_ratio.mat','amplitude_ratio');
% save('deviation_theta.mat','deviation_theta');
% 
% % load('result\厚度变化结果\边缘251,150\amplitude_ratio.mat');
% % load('result\厚度变化结果\边缘251,150\deviation_theta.mat');
% % load('amplitude_ratio.mat');
% % load('deviation_theta.mat');
% 
% figure
% plot(d,amplitude_ratio.*100,'LineWidth',2);
% grid on
% title('加罩前后点源成像幅值比/%');
% 
% figure
% plot(d,deviation_theta,'LineWidth',2);
% grid on
% title('加罩后点源成像偏离角/度');

% 介电常数影响
epsilon_relative = 2:0.1:5;
amplitude_ratio = zeros(1,length(epsilon_relative));
deviation_theta = zeros(1,length(epsilon_relative));
parfor kk = 1 : length(epsilon_relative)
    [~, ~, ~, ~,  amp_ratio, delta_theta ] = GetResult2(epsilon_relative(kk))
    amplitude_ratio(kk) = amp_ratio;
    deviation_theta(kk) = delta_theta;
end

save('amplitude_ratio.mat','amplitude_ratio');
save('deviation_theta.mat','deviation_theta');
% load('amplitude_ratio.mat');
% load('deviation_theta.mat');


figure
plot(epsilon_relative,amplitude_ratio.*100,'LineWidth',2);
grid on
xlabel('天线罩介质相对介电常数'); 
ylabel('加罩前后点源成像幅值比/%');
title('加罩前后点源成像幅值比随天线罩介质相对介电常数变化关系');

figure
plot(epsilon_relative,deviation_theta,'LineWidth',2);
grid on
xlabel('天线罩介质相对介电常数'); 
ylabel('加罩后点源成像偏离角/度');
title('加罩后点源成像偏离角随天线罩介质相对介电常数变化关系');


toc