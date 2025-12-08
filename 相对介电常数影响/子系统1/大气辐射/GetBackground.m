function [radio_atmosparam,target_atmosparam,targetdown_angle] = GetBackground( Param )
% 作用：计算无目标时背景的辐射亮温。
%输出观测点所在高度altitude
%输出不同仰角亮温分布的拟合系数fitting_coe
% fit_num = 100;      % 拟合点的数目
normal_vector = [0 0 1];        % 地面法向量 
polar_vector = Param.Observe.polar_vector;          % 极化向量
% hwait = waitbar(0,'参数初始化');             % 初始化进度条
%% 辐射计处背景亮温求解
elevation_angle = linspace(-90,90,1801);
% BRITEMP = elevation_angle;
radio_altitude = Param.Observe.radio_altitude;

radio_atmosparam = zeros(18,100);
parfor_progress(size(radio_atmosparam,1));                    % 用于记录并行计算进度
parfor m = 1:18
    pitch_angle = elevation_angle;
    BRITEMP = radio_atmosparam(m,:);
    for n = 1:100
        % 反推路径，计算亮温、衰减、出射角等
        [mode,BRITEMP(n),ATT,Beta] = GetAtmosRadiation( radio_altitude, pitch_angle((m-1)*100+n),Param );
        % 如果路径遇到大地，即mod==1
        if mode == 1
            [~,BT_re,~,~] = GetAtmosRadiation( 0, Beta-90,Param );  % Beta-90为路径在地面反射向量的俯仰角
            inc_vector = [-cosd(Beta-90),0,-sind(Beta-90)];          % 入射向量
            [eh_ground,ev_ground] = RefGround(180-Beta,Param.Ambient);            % 计算地面发射率分量
            e_ground = LinearPolarization(polar_vector,inc_vector,normal_vector,eh_ground,ev_ground);     % 计算地面发射率
            BT_grd = BT_re*(1-e_ground) + Param.Ambient.ground_temp*e_ground;
            BRITEMP(n) = BRITEMP(n) + BT_grd*exp(-0.23*ATT);
            %     else break;
        end
    end
    radio_atmosparam(m,:) = BRITEMP;
     parfor_progress;
end
parfor_progress(0);
radio_atmosparam = reshape(radio_atmosparam.',[1,18*100]);
%% 计算入射角为90度时辐射计处大气辐射
[mode,BRITEMP_90,ATT,Beta] = GetAtmosRadiation( radio_altitude, 90.0,Param );
% 如果路径遇到大地，即mod==1
if mode == 1
    [~,BT_re,~,~] = GetAtmosRadiation( 0, Beta-90,Param );  % Beta-90为路径在地面反射向量的俯仰角
    inc_vector = [-cosd(Beta-90),0,-sind(Beta-90)];          % 入射向量
    [eh_ground,ev_ground] = RefGround(180-Beta,Param.Ambient);            % 计算地面发射率分量
    e_ground = LinearPolarization(polar_vector,inc_vector,normal_vector,eh_ground,ev_ground);     % 计算地面发射率
    BT_grd = BT_re*(1-e_ground) + Param.Ambient.ground_temp*e_ground;
    BRITEMP_90 = BRITEMP_90 + BT_grd*exp(-0.23*ATT);
    %     else break;
end
radio_atmosparam = [radio_atmosparam,BRITEMP_90,0];

% radiodown_angle = elevation_angle(m-1);                 % 辐射计处向下观测时反射地面临界角
% elevation_angle = elevation_angle(1:m-1);               % 取向下反射地面的俯仰角
% BRITEMP = BRITEMP(1:m-1);
% radiodown_fitting_coeff = polyfit(elevation_angle,BRITEMP,10);      % 十次多项式拟合，radiodown_fitting_coeff为每项的拟合系数向量，包含11个元素

% figure,plot(elevation_angle,BRITEMP,elevation_angle,polyval(radiodown_fitting_coeff,elevation_angle))
% figure,plot(elevation_angle,BRITEMP)
% figure,plot(elevation_angle,polyval(radiodown_fitting_coeff,elevation_angle))
%% 目标处亮温向下拟合系数求解
target_altitude = Param.Observe.target_altitude;
elevation_angle = linspace(-90,90,1801);
target_atmosparam = zeros(3,1802);
BRITEMP = zeros(18,100);
ATT = zeros(18,100);
BT_re = zeros(18,100);
parfor_progress(size(BRITEMP,1));                    % 用于记录并行计算进度
parfor m = 1:18
    pitch_angle = elevation_angle;
    BRITEMP1 = BRITEMP(m,:);
    ATT1 = ATT(m,:);
    BT_re1 = BT_re(m,:);
    for n = 1:100
        % 反推路径，计算亮温、衰减、出射角等
        [mode,BRITEMP1(n),ATT1(n),Beta] = GetAtmosRadiation( target_altitude, pitch_angle((m-1)*100+n),Param );
        
        % 如果路径遇到大地，即mod==1
        if mode == 1
            [~,BT_re1(n),~,~] = GetAtmosRadiation( 0, Beta-90,Param );  % Beta-90即对大地的入射角
            
        end
    end
    BRITEMP(m,:) = BRITEMP1;
    ATT(m,:) = ATT1;
    BT_re(m,:) = BT_re1;
    parfor_progress;
end
parfor_progress(0);

%% 计算入射角为90度时目标处大气辐射

% 反推路径，计算亮温、衰减、出射角等
[mode,BRITEMP_90,ATT_90,Beta] = GetAtmosRadiation( target_altitude, 90.0,Param );
BT_re_90 = 0;
% 如果路径遇到大地，即mod==1
if mode == 1
    [~,BT_re_90,~,~] = GetAtmosRadiation( 0, Beta-90,Param );  % Beta-90即对大地的入射角
    
end

target_atmosparam(1,:) = [reshape(BRITEMP.',[1,18*100]),BRITEMP_90,0];
target_atmosparam(2,:) = [reshape(ATT.',[1,18*100]),ATT_90,0];
target_atmosparam(3,:) = [reshape(BT_re.',[1,18*100]),BT_re_90,0];
count = find(target_atmosparam(3,:),1,'last');
targetdown_angle = elevation_angle(count);                % 目标处向下观测时反射地面临界角
target_atmosparam(2,count+1) = 0;

% figure,plot(elevation_angle,BRITEMP,elevation_angle,polyval(targetdown_fitting_coeff(1,:),elevation_angle))
% figure,plot(elevation_angle,ATT,elevation_angle,polyval(targetdown_fitting_coeff(2,:),elevation_angle))
% figure,plot(elevation_angle,BT_re,elevation_angle,polyval(targetdown_fitting_coeff(3,:),elevation_angle))
% figure,plot(elevation_angle,polyval(targetdown_fitting_coeff(1,:),elevation_angle))
% figure,plot(elevation_angle,BRITEMP)
% figure,plot(elevation_angle,polyval(targetdown_fitting_coeff(2,:),elevation_angle))
% figure,plot(elevation_angle,ATT)
% figure,plot(elevation_angle,polyval(targetdown_fitting_coeff(3,:),elevation_angle))
% figure,plot(elevation_angle,BT_re)
%% 辐射计处中段拟合系数求解
% elevation_angle = linspace(radiodown_angle,radiodown_angle+2,floor(fit_num/4));
% BRITEMP = elevation_angle;
% for m = 1:length(elevation_angle)
%     % 反推路径，计算亮温、衰减、出射角等
%     [~,BRITEMP(m),~,~] = GetAtmosRadiation( Param.Observe.target_altitude, elevation_angle(m),Param );
%     
%     % 如果路径遇到大地，即mod==1
% %     if mode == 1
% %         [~,BT_re(m),~,~] = GetAtmosRadiation( 0, Beta-90,Param );  % Beta-90即对大地的入射角
% %         inc_vector = [-cosd(Beta-90),0,-sind(Beta-90)];          % 入射向量 
% %         [eh_ground,ev_ground] = RefGround(180-Beta,Param.Ambient);            % 计算地面反射率分量
% %         [e_ground,~] = LinearPolarization(polar_vector,inc_vector,normal_vector,eh_ground,ev_ground);     % 计算地面发射率
% %         BT_grd = BT_re(m)*(1-e_ground) + Param.Ambient.ground_temp*e_ground;
% %         BRITEMP(m) = BRITEMP(m) + BT_grd*exp(-0.23*ATT(m));
% %     end
%     waitbar((m+4*fit_num)/fit_num/6.5,hwait,'参数初始化');
% end
% radiomid_fitting_coeff = polyfit(elevation_angle,BRITEMP,15);      % 十五次多项式拟合，radiomid_fitting_coeff为每项的拟合系数向量，包含16个元素

% figure,plot(elevation_angle,BRITEMP,elevation_angle,polyval(radiomid_fitting_coeff,elevation_angle))
% figure,plot(elevation_angle,BRITEMP)
% figure,plot(elevation_angle,polyval(radiomid_fitting_coeff,elevation_angle))
%% 目标处中段拟合系数求解
% elevation_angle = linspace(targetdown_angle,targetdown_angle+2,floor(fit_num/2));
% BRITEMP = elevation_angle;
% for m = 1:length(elevation_angle)
%     % 反推路径，计算亮温、衰减、出射角等
%     [~,BRITEMP(m),~,~] = GetAtmosRadiation( Param.Observe.target_altitude, elevation_angle(m),Param );
%     
%     % 如果路径遇到大地，即mod==1
% %     if mode == 1
% %         [~,BT_re(m),~,~] = GetAtmosRadiation( 0, Beta-90,Param );  % Beta-90即对大地的入射角
% %         inc_vector = [-cosd(Beta-90),0,-sind(Beta-90)];          % 入射向量 
% %         [eh_ground,ev_ground] = RefGround(180-Beta,Param.Ambient);            % 计算地面反射率分量
% %         [e_ground,~] = LinearPolarization(polar_vector,inc_vector,normal_vector,eh_ground,ev_ground);     % 计算地面发射率
% %         BT_grd = BT_re(m)*(1-e_ground) + Param.Ambient.ground_temp*e_ground;
% %         BRITEMP(m) = BRITEMP(m) + BT_grd*exp(-0.23*ATT(m));
% %     end
%     waitbar((m+2000)/2150,hwait,'参数初始化');
% end                
% % 11次多项式拟合，targetmid_fitting_coeff为每项的拟合系数向量，包含12个元素,分别对BRITEMP、ATT、BT_re进行拟合
% targetmid_fitting_coeff = polyfit(elevation_angle,BRITEMP,15);        % targetmid_fitting_coeff(1,:)对BRITEMP进行拟合

% figure,plot(elevation_angle,BRITEMP,elevation_angle,polyval(targetmid_fitting_coeff,elevation_angle))
% figure,plot(elevation_angle,polyval(targetmid_fitting_coeff,elevation_angle))
% figure,plot(elevation_angle,BRITEMP)
%% 辐射计处向上拟合系数求解
% elevation_angle = linspace(radiodown_angle+2,90,fit_num);
% BRITEMP = elevation_angle;
% for m = 1:fit_num
%     % 反推路径，计算亮温、衰减、出射角等
%     [~,BRITEMP(m),~,~] = GetAtmosRadiation( Param.Observe.radio_altitude, elevation_angle(m),Param );
%     waitbar((m+4.5*fit_num)/fit_num/6.5,hwait,'参数初始化');
% end
% radioup_fitting_coeff = polyfit(elevation_angle,BRITEMP,15);      % 十次多项式拟合，radioup_fitting_coeff为每项的拟合系数向量，包含11个元素

% figure,plot(elevation_angle,BRITEMP,elevation_angle,polyval(radioup_fitting_coeff,elevation_angle))
%  figure,plot(elevation_angle,BRITEMP)
%  figure,plot(elevation_angle,polyval(radioup_fitting_coeff,elevation_angle))
 %% 目标处向上拟合系数求解
%  elevation_angle = linspace(targetdown_angle+2,90,200);
%  BRITEMP = elevation_angle;
% for m = 1:200
%     % 反推路径，计算亮温、衰减、出射角等
%     [~,BRITEMP(m),~,~] = GetAtmosRadiation( Param.Observe.target_altitude, elevation_angle(m),Param );
% %     waitbar((m+2050)/2150,hwait,'参数初始化');
% end
% targetup_fitting_coeff = polyfit(elevation_angle,BRITEMP,15);      %十次多项式拟合，targetup_fitting_coeff为每项的拟合系数向量，包含11个元素
% 
% figure,plot(elevation_angle,BRITEMP,elevation_angle,polyval(targetup_fitting_coeff,elevation_angle))
% figure,plot(elevation_angle,polyval(targetup_fitting_coeff,elevation_angle))
% figure,plot(elevation_angle,BRITEMP)
% pause(0.0001)
% close(hwait)
end

