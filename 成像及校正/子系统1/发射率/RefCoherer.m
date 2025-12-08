function [eh,ev] = RefCoherer(incidence,Param)
%计算金属板的反射率
%输入：探测方向矢量与面方向矢量的夹角incidence
permittivity =  Param.Ambient.atmos_permittivity;       % 自由空间介电常数
frequency = Param.Observe.center_freq;                  % 辐射计中心频率
conductivity = Param.Target.target_conductivity;        % 目标电导率

%入射角incidence为辐射射线金属平面方向矢量的夹角
t=sqrt(4.*pi.*frequency*1e9.*permittivity./conductivity);
Rh=1-2.*t.*cosd(incidence); %反射率水平极化分量
Rv=(2.*(cosd(incidence).^2)-2.*t.*cosd(incidence)+t.^2)./(2.*(cosd(incidence).^2)+2.*t.*cosd(incidence)+t.^2);
%反射率垂直极化分量
eh=1-Rh;        % 由反射分量得到发射分量
ev=1-Rv;
end