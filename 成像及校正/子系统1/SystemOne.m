function [ TB ] = SystemOne( tt )
%SYSTEMONE 此处显示有关此函数的摘要
%   此处显示详细说明
global SimuParam;

SetParamOne(tt);
    
[TB,~,~,~,~] = AccGetRadiationTemp(SimuParam);           % 计算亮温分布
    
colorbar;
figure,imagesc(TB)
axis image;

end

