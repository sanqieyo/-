function [mode,BRITEMP,ATT,Beta] = GetAtmosRadiation( altitude, elevation_angle,Param)
%该函数的作用：从目标位置，根据仰角反追踪辐射路径，当路径穿出最外层，或遇到大地，停止追踪

%输入观测高度altitude，观测仰角elevation_angle，结构体Param
%输出结束模式mode，大气辐射BRITEMP，大气衰减ATT，出射角Beta


% 相关数据载入：
R = Param.Ambient.earth_radius ; 
Top = Param.Atmos.Top;
Thickness = Param.Atmos.Thickness;
BriTemp = Param.Atmos.BriTemp;
Att = Param.Atmos.Att;
if elevation_angle==0
end

% 计算观测点所在层
Layer = 0; % 记录层数
for i = 1:Param.Atmos.LayerNum
    if Top(i)-R > altitude
        Layer = i;
        break;
    end
end

% 以下是计算过程
Alpha = 90 - elevation_angle; %记录每一层的入射角
HO = altitude + R; % 记录入射点离地球球心的距离
angle = 0; % 记录路径的偏心角
ATT = 0; % 记录路径的衰减因子
BRITEMP = 0; % 记录各层大气辐射的贡献
length = 0; % 记录路径长度
while 1
    if Alpha == 0       % 完全向上看
        Beta = 0;       % 记录每一层的出射角
        Distance = Top(Layer) - HO;
        length = length + Distance;
        BRITEMP = BRITEMP + ...
            BriTemp(Layer)*...
            exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
        ATT = ATT + Att(Layer)*Distance;
        if Layer == Param.Atmos.LayerNum    % 路径到达最外层
            HOO = Top(Param.Atmos.LayerNum);        % 用于记录终止点（目标、大地、最外层）离地球球心的距离
            mode = 2;           % 路径终止模式：1.遇到大地 / 2.穿过最外层 / 3.遇到目标
            break;
        end
        HO = Top(Layer);
        Layer = Layer + 1;
    elseif Alpha == 180  % 完全向下看
        Beta = 180;         % 记录每一层的出射角
        Distance = -Top(Layer)+Thickness(Layer) + HO;
        length = length + Distance;
        BRITEMP = BRITEMP + ...
            BriTemp(Layer)*...
            exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
        ATT = ATT + Att(Layer)*Distance;
        if Layer == 1  % 路径到达地面
            HOO = R;
            mode = 1;
            break;
        end
        HO = Top(Layer-1);
        Layer = Layer - 1;
    elseif 180-asind((Top(Layer)-Thickness(Layer))/HO)>=Alpha  % 情况一：路径向外层拓展
        Beta = asind(sind(Alpha)*HO/Top(Layer));
        Gamma = Alpha - Beta;
        angle = angle + Gamma;
        Distance = Top(Layer)*sind(Gamma)/sind(Alpha);
        length = length + Distance;
        BRITEMP = BRITEMP + ...
            BriTemp(Layer)*...
            exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
        ATT = ATT + Att(Layer)*Distance;
        if Layer == Param.Atmos.LayerNum  % 路径到达最外层
            HOO = Top(Param.Atmos.LayerNum);
            mode = 2;
            break;
        end
        Alpha = asind(sind(Beta)*Param.Atmos.Refr(Layer)/Param.Atmos.Refr(Layer+1));
        if imag(Alpha) ~= 0  % 路径遇到全反射陷阱，以下进行特殊化处理
%             warning('仿真程序遇到全反射情况，可能带来较大仿真误差！！！');
            Alpha = 90;
        end
        HO = Top(Layer);
        Layer = Layer + 1;
    else  %情况二：路径向内层拓展
        Beta = 180 - asind(sind(Alpha)*HO/(Top(Layer)-Thickness(Layer)));
        Gamma = Alpha - Beta;
        angle = angle + Gamma;
        Distance = (Top(Layer)-Thickness(Layer))*sind(Gamma)/sind(Alpha);
        length = length + Distance;
        BRITEMP = BRITEMP + ...
            BriTemp(Layer)*...
            exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
        ATT = ATT + Att(Layer)*Distance;
        if Layer == 1  % 路径到达地面
            HOO = R;
            mode = 1;
            break;
        end
        Alpha = 180 - ...
            asind(sind(Beta).*Param.Atmos.Refr(Layer)./Param.Atmos.Refr(Layer-1));
        HO = Top(Layer-1);
        Layer = Layer - 1;
    end
end
end

