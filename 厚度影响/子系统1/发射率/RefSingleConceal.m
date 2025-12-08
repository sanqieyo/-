function [eh,ev] = RefSingleConceal(incidence,Param)
%计算隐身目标的发射率
%输入：入射角incidence
%输出：输出发射率的水平极化eh和垂直极化分量ev

Frequency = Param.Observe.center_freq;                 % 辐射计中心频率
Permittivity =  Param.Target.layer1_target_permittivity;      % 目标吸波材料的介电常数
Permeability = Param.Target.layer1_target_permeability;       % 目标吸波材料的磁导率
Thick = Param.Target.layer1_coat_thickness;                   % 目标吸波材料的厚度，单位m
Wavelength = 0.3./Frequency;

ee = exp((-1i.*4.*Thick.*pi.*Permittivity.*Permeability)./(Wavelength.*sqrt(Permittivity.*Permeability-(sind(incidence)).^2)));
r1 = Permeability.*cosd(incidence) - sqrt(Permeability.*Permittivity - (sind(incidence)).^2);
r11 = Permeability.*cosd(incidence) + sqrt(Permeability.*Permittivity - (sind(incidence)).^2); 
rr = r1./r11;       % 反射系数R1
Rh = abs((rr - ee)./(1 - rr.*ee)).^2;

eh=1-Rh;

ee = exp((-1i.*4.*Thick.*pi.*Permittivity.*Permeability)./(Wavelength.*sqrt(Permittivity.*Permeability-(sind(incidence)).^2)));
r1 = Permittivity.*cosd(incidence) - sqrt(Permeability.*Permittivity - (sind(incidence)).^2);
r11 = Permittivity.*cosd(incidence) + sqrt(Permeability.*Permittivity - (sind(incidence)).^2); 
rr = r1./r11;
Rv = abs((rr + ee)./(1+rr.*ee)).^2;
ev=1-Rv;

end