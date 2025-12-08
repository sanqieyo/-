function Param = ProcessAtmosParam(Param)

% 该函数用于设定模型参数：
%     1.无参数时，对模型参数进行默认初始化
%     2.有参数（两个）时，改变参数设定
%       例：若想改变SimuParam.Radio.WorkFreq为118，则调用形式为：SetParam('Radio.WorkFreq',118)

workfreq = Param.Observe.center_freq;
Param.Atmos.LayerNum = 918;
Param.Atmos.Thickness = logspace(-4, 0, Param.Atmos.LayerNum);  % 大气每层的厚度（km）
Param.Atmos.Top = zeros(1, Param.Atmos.LayerNum);               % 大气各层顶部到地球球心的距离（km）
Param.Atmos.Temp = zeros(1, Param.Atmos.LayerNum);              % 大气每层物理温度（K）
Param.Atmos.Ke = zeros(1, Param.Atmos.LayerNum);                % 大气消光系数（db/km）
Param.Atmos.Att = zeros(1, Param.Atmos.LayerNum);               % 大气衰减因子（db/km）
Param.Atmos.Emi = zeros(1, Param.Atmos.LayerNum);               % 大气发射系数（db/km）
Param.Atmos.Refr = zeros(1, Param.Atmos.LayerNum);              % 大气折射率
Param.Atmos.BriTemp = zeros(1, Param.Atmos.LayerNum);           % 大气辐射亮温（db/km）
bottom = 0;
for m = 1 : Param.Atmos.LayerNum
    [P,T,ro] = StdAtmos(Param.Ambient.atmos_type, bottom+Param.Atmos.Thickness(m)/2,Param.Ambient.ground_temp);
    Param.Atmos.Temp(m) = T;
    Param.Atmos.Ke(m) = itu676a(workfreq, P, T, ro) + GetWeatherAtt(Param.Ambient.weather_type, T, ...
        workfreq,bottom+Param.Atmos.Thickness(m)/2);
    Param.Atmos.Att(m) = Param.Atmos.Ke(m);
    Param.Atmos.Emi(m) = 1-10^(-Param.Atmos.Ke(m)./10);         %发射系数
    Param.Atmos.BriTemp(m) = Param.Atmos.Temp(m).*Param.Atmos.Emi(m);   %每Km大气的辐射亮温
    N = RefAtmos(P,T,ro)+itu676d(workfreq,P,T,ro);
    Param.Atmos.Refr(m) = 1+N.*(1e-6);
    bottom = bottom + Param.Atmos.Thickness(m);
    Param.Atmos.Top(m) = bottom + Param.Ambient.earth_radius;
end
end




