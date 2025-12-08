function [ lambda, Gain, Omega] = ProcessSystemInput( Param )
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
lambda = Param.c / Param.f0;% 波长
Gain = 10^(Param.gain / 10);
Omega = 4*pi / Gain;

end

