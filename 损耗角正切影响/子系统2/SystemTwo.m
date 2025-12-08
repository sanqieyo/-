function [ Visib, visibility_FF ,visibility_NF, visibility_Ideal] = SystemTwo( TB )
% 子系统2运行函数 
% 输入：    TB                   子系统1生成的亮温矩阵

% 输出：    Visib                头罩影响下的可见度

%           visibility_FF        头罩幅相误差影响的可见度
%           visibility_NF        头罩热辐射影响可见度
%           visibility_Ideal     没有头罩影响的理想可见度

addpath(genpath(pwd))
global SimuParam;
SetParamTwo();

[visibility_FF] = GetFFVisibility( SimuParam ,TB);
[visibility_NF] = GetNFVisibility( SimuParam );
Visib = visibility_FF + visibility_NF;    % 无冗余可见度输出

[visibility_Ideal] = GetIdealVisibility( SimuParam ,TB);

save('visibility_FF','visibility_FF');
save('visibility_NF','visibility_NF');
save('Visib','Visib');
save('visibility_Ideal','visibility_Ideal');

end

