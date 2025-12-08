function [BRITEMP,ATT,mode] = GetAtmosAttenuation( Param )
%该函数的作用：从辐射计位置，根据仰角反追踪辐射路径，当路径遇到目标停止mode = 1，或遇到最外层、大地报错停止mode = 0。

%输入观测高度altitude，观测仰角elevation_angle，结构体Param
%输出结束模式mode，大气辐射BRITEMP，大气衰减ATT
%输出终止点（目标、大地、最外层）离地球球心的距离HOO，出射角Beta，路径长度length，路径偏心角angle
% 相关数据载入：
altitude = Param.Observe.radio_altitude;                %辐射计所在高度
elevation_angle = Param.Observe.elevation_angle;        %观测俯仰角
distance = Param.Observe.observe_dist;                  %观测距离
R = Param.Ambient.earth_radius ;                        %地球半径
Top = Param.Atmos.Top;
Thickness = Param.Atmos.Thickness;
BriTemp = Param.Atmos.BriTemp;
Att = Param.Atmos.Att;
% 计算观测点所在层
Layer = 0;      % 记录层数
for i = 1:Param.Atmos.LayerNum
    if Top(i)-R > altitude
        Layer = i;
        break;
    end
end

% 以下是计算过程
Beta = 0;       % 记录每一层的
Alpha = 90 - elevation_angle;       %记录每一层的入射角
HO = altitude + R;                  % 记录入射点离地球球心的距离
angle = 0;      % 记录路径的偏心角
ATT = 0;        % 记录路径的衰减因子
BRITEMP = 0;    % 记录各层大气辐射的贡献
length = 0;     % 记录路径长度
HOO = 0;        % 用于记录终止点（目标、大地、最外层）离地球球心的距离
while 1
    if Alpha == 0  % 完全向上看
        Beta = 0;
        Distance = Top(Layer) - HO;
        length = length + Distance;
        if length >= distance       %路径到达目标
            Distance = Distance - ( length - distance );
            HOO = R + altitude + distance;      %用于记录终止点（目标、大地、最外层）离地球球心的距离
            BRITEMP = BRITEMP + ...
                BriTemp(Layer)*...
                exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
            ATT = ATT + Att(Layer)*Distance;
            length = distance;
            mode = 1;
            break;
        end
        BRITEMP = BRITEMP + ...
            BriTemp(Layer)*...
            exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
        ATT = ATT + Att(Layer)*Distance;
        if Layer == Param.Atmos.LayerNum  % 路径到达最外层
            mode = 0;
            HOO = Top(Param.Atmos.LayerNum);
            break;
        end
        HO = Top(Layer);
        Layer = Layer + 1;
    elseif Alpha == 180         % 完全向下看
        Beta = 180;
        Distance = -Top(Layer)+Thickness(Layer) + HO;
        length = length + Distance;
        if length >= distance   %路径到达目标
            Distance = Distance - ( length - distance );
            HOO = R + altitude - distance;
            BRITEMP = BRITEMP + ...
                BriTemp(Layer)*...
                exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
            ATT = ATT + Att(Layer)*Distance;
            length = distance;
            mode = 1;
            break;
        end
        BRITEMP = BRITEMP + ...
            BriTemp(Layer)*...
            exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
        ATT = ATT + Att(Layer)*Distance;
        if Layer == 1           % 路径到达地面
            if abs( length-distance ) < 1e-4
                HOO = R;
                mode = 1;
                break;
            end
            mode = 0;
            HOO = R;
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
        if length >= distance  %路径到达目标
            Distance = Distance - ( length - distance );
            HOO = ( HO^2+Distance^2+2*HO*Distance*abs(cosd(Alpha)) )^0.5;
            Beta = asind( sind(Alpha)*HO/HOO );
            Gamma = asind( sind(Alpha)*Distance/HOO );
            angle = angle + Gamma;
            BRITEMP = BRITEMP + ...
                BriTemp(Layer)*...
                exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
            ATT = ATT + Att(Layer)*Distance;
            length = distance;
            mode = 1;
            break;
        end
        BRITEMP = BRITEMP + ...
            BriTemp(Layer)*...
            exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
        ATT = ATT + Att(Layer)*Distance;
        if Layer == Param.Atmos.LayerNum  % 路径到达最外层
            mode=4;
            HOO = Top(Param.Atmos.LayerNum);
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
        if length >= distance
            Distance = Distance - ( length - distance );
            HOO = ( HO^2+Distance^2-2*HO*Distance*abs(cosd(Alpha)) )^0.5;
            Beta = 180 - asind( sind(Alpha)*HO/HOO );
            Gamma = asind( sind(Alpha)*Distance/HOO );
            angle = angle + Gamma;
            BRITEMP = BRITEMP + ...
                BriTemp(Layer)*...
                exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
            ATT = ATT + Att(Layer)*Distance;
            length = distance;
            mode = 1;
            break;
        end
        BRITEMP = BRITEMP + ...
            BriTemp(Layer)*...
            exp(-0.23*ATT)/0.23/Att(Layer)*(1-exp(-0.23*Att(Layer)*Distance));
        ATT = ATT + Att(Layer)*Distance;
        if Layer == 1  % 路径到达地面
            if abs( length-distance ) < 1e-4
                HOO = R;
                mode = 1;
                break;
            end
            mode = 0;
            HOO = R;
            break;
        end
        Alpha = 180 - ...
            asind(sind(Beta).*Param.Atmos.Refr(Layer)./Param.Atmos.Refr(Layer-1));
        HO = Top(Layer-1);
        Layer = Layer - 1;
    end
end
end

