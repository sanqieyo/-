function [TB,e_tg,count,target_area,TB_target] = AccGetRadiationTemp(Param)
%本函数用于计算复杂场景的亮温分布

y_list = Param.Observe.horizon_scale;           % 横向扫描序列
z_list = Param.Observe.vertical_scale;          % 纵向扫描序列

target_point = Param.Target.boat.target_point;       % 原始坐标系中的目标节点坐标
trans_point = target_point*Param.Observe.rotate_factor;     % 旋转后坐标系中的目标节点坐标
zenith_vector = [0;0;1];        % 铅垂方向向量

TB = zeros(length(z_list),length(y_list));      % 亮温矩阵初始化
e_tg = TB;
count = TB;

radio_atmos = Param.RadTemp.radio_atmosparam;
%% 背景和目标分开计算
if Param.calculation_mode == 0                  % 常规计算
    observe_dist = Param.Observe.observe_dist*1000;
    
    % 计算包围盒范围
    x_node = [max(trans_point(:,1)),min(trans_point(:,1))]-observe_dist;
    y_node = -[max(trans_point(:,2)),min(trans_point(:,2))]*observe_dist;
    z_node = -[max(trans_point(:,3)),min(trans_point(:,3))]*observe_dist;
    box_y = y_node'*(1./x_node);
    box_z = z_node'*(1./x_node);
    z_max = max(max(box_z));
    y_max = max(max(box_y));
    z_min = min(min(box_z));
    y_min = min(min(box_y));
    
    % 目标区域亮温计算
    hwait = waitbar(0,'目标区域亮温计算');
    waitbar_length = ceil((atand((y_max)/observe_dist)-atand((y_min)/observe_dist))/Param.Observe.pixel_space);
    % 目标区域亮温计算进度条长度
    waitbar_zero = length(y_list)/2+floor(atand((y_min)/observe_dist)/Param.Observe.pixel_space);         % 目标区域亮温计算进度条起始点
    for m=1:length(y_list)
        y_loc = observe_dist*tand(y_list(m));                   % 像素点在成像平面上的y坐标
        if y_loc<=y_max&&y_loc>=y_min                       % 确定目标区域y坐标大致范围
            for n = 1:length(z_list)
                z_loc = tand(z_list(n))*observe_dist/cosd(y_list(m));      % 像素点在成像平面上的z坐标
                if z_loc<=z_max&&z_loc>=z_min               % 确定目标区域y坐标大致范围
                    scan_vector = [-observe_dist, y_loc , z_loc];            % 旋转后的坐标系中入射射线的方向矢量
                    inc_vector = scan_vector/Param.Observe.rotate_factor;           % 入射射线在原始坐标系的方向矢量
                    a = sind(z_list(n));
                    b = cosd(z_list(n));
                    c = sind(y_list(m));
                    d = cosd(y_list(m));
                    rotate_factor_polar=[b*d c -a*d
                        -b*c d a*c
                        a  0  b];                    % 坐标转换矩阵，由scan_vector旋转到[-observe_dist, 0 , 0]
                    polar_vector_scan = Param.Observe.polar_vector/rotate_factor_polar/Param.Observe.rotate_factor;         % 不同观测方向上的极化矢量
                    [TB(n,m),count(n,m),e_tg(n,m)] = FindFacet(inc_vector,Param.Observe.observe_position,Param,polar_vector_scan);           % 计算像素点对应目标面元的亮温值
                    if TB(n,m)==0               % 若像素点不在目标上，即入射射线与目标不相交，直接应用拟合系数计算面元亮温
                        zenith_theta = acosd(inc_vector*zenith_vector/norm(inc_vector));                % 计算入射射线与地表法向量的夹角，即为天顶角
                        zenith_theta = 90 - zenith_theta;
                        index = ceil((zenith_theta-(-90)+0.00001)/0.1);              % 判断俯仰角所在区间
                        factor = ceil((zenith_theta-(-90)+0.00001)/0.1)-(zenith_theta-(-90)+0.00001)/0.1;                        % 加权因子
                        TB(n,m) = radio_atmos(index)*factor+radio_atmos(index+1)*(1-factor);        % 加权计算对应俯仰角下的背景辐射
                    end
                end
            end
            pause(0.1)
            waitbar((m-waitbar_zero)/waitbar_length,hwait,'目标区域亮温计算');
        end
    end
    
    % 背景区域亮温计算
    for m = 1:length(z_list)
        for n = 1:length(y_list)
            y_loc = observe_dist*tand(y_list(n));                   % 像素点在成像平面上的y坐标
            z_loc = tand(z_list(m))*observe_dist/cosd(y_list(n));      % 像素点在成像平面上的z坐标
            if TB(m,n)==0       %若像素点不在目标区域内
                scan_vector = [-observe_dist, y_loc , z_loc];        % 旋转后的坐标系中入射射线的方向矢量
                inc_vector = scan_vector/Param.Observe.rotate_factor;           % 入射射线在原始坐标系的方向矢量
                zenith_theta = acosd(inc_vector*zenith_vector/norm(inc_vector));                % 计算入射射线与地表法向量的夹角，即为天顶角
                zenith_theta = 90 - zenith_theta;
                ind ex = ceil((zenith_theta-(-90)+0.00001)/0.1);              % 判断俯仰角所在区间
                factor = ceil((zenith_theta-(-90)+0.00001)/0.1)-(zenith_theta-(-90)+0.00001)/0.1;                        % 加权因子
                TB(m,n) = radio_atmos(index)*factor+radio_atmos(index+1)*(1-factor);        % 加权计算对应俯仰角下的背景辐射
            end
        end
        waitbar(m/length(z_list),hwait,'背景区域亮温计算');
    end
    close(hwait)
    pixel_area = (observe_dist*tand(Param.Observe.pixel_space))^2;
    target_area = pixel_area*nnz(count);
    TB_target = sum(TB(count~=0))/nnz(count);                        % 目标区域平均亮温
end

%% 并行计算（背景和目标一起计算）
if  Param.calculation_mode == 1                  % 并行计算
    observe_dist = Param.Observe.observe_dist*1000;
    rotate_factor = Param.Observe.rotate_factor;
%     polar_angle = Param.Observe.polar_angle;    % 极化角

    % 计算包围盒范围
    x_node = [max(trans_point(:,1)),min(trans_point(:,1))]-observe_dist;
    y_node = -[max(trans_point(:,2)),min(trans_point(:,2))]*observe_dist;
    z_node = -[max(trans_point(:,3)),min(trans_point(:,3))]*observe_dist;
    box_y = y_node'*(1./x_node);
    box_z = z_node'*(1./x_node);
    z_max = max(max(box_z));
    y_max = max(max(box_y));
    z_min = min(min(box_z));
    y_min = min(min(box_y));

    observe_position = Param.Observe.observe_position;
    polar_vector = Param.Observe.polar_vector;
    y_length = length(y_list);
    parfor_progress(length(z_list));                    % 用于记录并行计算进度

    parfor m=1:length(z_list)
        radio_atmos1 = radio_atmos;
        TB1 =zeros(1,y_length);
        e_tg1=TB1;
        count1=TB1;
        for n = 1:length(y_list)
            y_loc = observe_dist*tand(y_list(n));                   % 像素点在成像平面上的y坐标
            z_loc = tand(z_list(m))*observe_dist/cosd(y_list(n));      % 像素点在成像平面上的z坐标
            scan_vector = [-observe_dist, y_loc , z_loc];        % 旋转后的坐标系中扫描射线的方向矢量
            inc_vector = scan_vector/rotate_factor;       % 扫描射线在原始坐标系的方向矢量

            a = sind(z_list(m));
            b = cosd(z_list(m));
            c = sind(y_list(n));
            d = cosd(y_list(n));
            rotate_factor_polar=[b*d c -a*d
                -b*c d a*c
                a  0  b];                    % 坐标转换矩阵，由scan_vector旋转到[-observe_dist, 0 , 0]
            polar_vector_scan = polar_vector/rotate_factor_polar/rotate_factor;         % 不同观测方向上的极化矢量

%             ph = cross(inc_vector,zenith_vector)/norm(cross(inc_vector,zenith_vector));
%             pv = cross(ph,inc_vector)/norm(cross(ph,inc_vector));
%             polar_vector_scan = ph*sind(polar_angle)+pv*cosd(polar_angle);    

            if y_loc<=y_max&&y_loc>=y_min&&z_loc<=z_max&&z_loc>=z_min
                
                [TB1(1,n),count1(1,n),e_tg1(1,n)] = FindFacet(inc_vector,observe_position,Param,polar_vector_scan);       % 计算目标船只在该点的亮温值
                %面元预筛选，筛选入射点周围的面元参与具体的求交算法
%                 [ screening_target_rank, screening_normal_vector, screening_rank_point ] = PreScreening( inc_vector, observe_position, Param);
%                 [TB1(1,n),count1(1,n),e_tg1(1,n)] = AccFindFacet(inc_vector,observe_position,Param,polar_vector_scan, screening_target_rank, screening_normal_vector, screening_rank_point);       % 计算像素点对应目标面元的亮温值
            end
            if TB1(1,n)==0               % 若像素点不在目标上，即入射射线与目标不相交，直接应用拟合系数计算面元亮温
%                 zenith_theta = acosd(inc_vector*zenith_vector/norm(inc_vector));                % 计算入射射线与地表法向量的夹角，即为天顶角
%                 zenith_theta = 90 - zenith_theta;
%                 index = ceil((zenith_theta-(-90)+0.00001)/0.1);              % 判断俯仰角所在区间
%                 factor = ceil((zenith_theta-(-90)+0.00001)/0.1)-(zenith_theta-(-90)+0.00001)/0.1;                        % 加权因子
%                 TB1(1,n) = radio_atmos1(index)*factor+radio_atmos1(index+1)*(1-factor);        % 加权计算对应俯仰角下的背景辐射
                [ screening_target_rank, screening_normal_vector, screening_rank_point ] = PreScreening( inc_vector, observe_position, Param);
                [TB1(1,n),~,~] = AccFindFacet(inc_vector,observe_position,Param,polar_vector_scan, screening_target_rank, screening_normal_vector, screening_rank_point);       % 计算像素点对应目标面元的亮温值
            end
            if TB1(1,n) == 0
                zenith_theta = acosd(inc_vector*zenith_vector/norm(inc_vector));                % 计算入射射线与地表法向量的夹角，即为天顶角
                zenith_theta = 90 - zenith_theta;
                index = ceil((zenith_theta-(-90)+0.00001)/0.1);              % 判断俯仰角所在区间
                factor = ceil((zenith_theta-(-90)+0.00001)/0.1)-(zenith_theta-(-90)+0.00001)/0.1;                        % 加权因子
                TB1(1,n) = radio_atmos1(index)*factor+radio_atmos1(index+1)*(1-factor);        % 加权计算对应俯仰角下的背景辐射
            end

        end
        TB(m,:) = TB1;
        e_tg(m,:) = e_tg1;
        count(m,:) = count1;
        parfor_progress;
    end
    parfor_progress(0);
    pixel_area = (observe_dist*tand(Param.Observe.pixel_space))^2;              % 单个像素点的投影面积
    target_area = pixel_area*nnz(count);                        % 目标区域投影面积
    TB_target = sum(TB(count~=0))/nnz(count);                        % 目标区域平均亮温
end

%% 并行计算（背景和目标分开计算）
if Param.calculation_mode == 2                  % 并行计算
    observe_dist = Param.Observe.observe_dist*1000;
    rotate_factor = Param.Observe.rotate_factor;
    
    % 计算包围盒范围
    x_node = [max(trans_point(:,1)),min(trans_point(:,1))]-observe_dist;
    y_node = -[max(trans_point(:,2)),min(trans_point(:,2))]*observe_dist;
    z_node = -[max(trans_point(:,3)),min(trans_point(:,3))]*observe_dist;
    box_y = y_node'*(1./x_node);
    box_z = z_node'*(1./x_node);
    z_max = max(max(box_z));
    y_max = max(max(box_y));
    z_min = min(min(box_z));
    y_min = min(min(box_y));
    observe_position = Param.Observe.observe_position;
    polar_vector = Param.Observe.polar_vector;
    
    waitbar_length = ceil((atand((y_max)/observe_dist)-atand((y_min)/observe_dist))/Param.Observe.pixel_space);
    parfor_progress(waitbar_length);
    parfor m=1:length(y_list)
        radio_atmos1 = radio_atmos;
        TB1 =zeros(length(z_list),1);
        e_tg1=TB1;
        count1=TB1;
        y_loc = observe_dist*tand(y_list(m));                   % 像素点在成像平面上的y坐标
        if y_loc<=y_max&&y_loc>=y_min                      % 确定目标区域z坐标大致范围
            for n = 1:length(z_list)
                z_loc = tand(z_list(n))*observe_dist/cosd(y_list(m));      % 像素点在成像平面上的z坐标
                if z_loc<=z_max&&z_loc>=z_min               % 确定目标区域y坐标大致范围
                    scan_vector = [-observe_dist, y_loc , z_loc];            % 旋转后的坐标系中入射射线的方向矢量
                    inc_vector = scan_vector/rotate_factor;           % 入射射线在原始坐标系的方向矢量
                    
                    a = sind(z_list(n));
                    b = cosd(z_list(n));
                    c = sind(y_list(m));
                    d = cosd(y_list(m));
                    rotate_factor_polar=[b*d c -a*d
                        -b*c d a*c
                        a  0  b];                    % 坐标转换矩阵，由scan_vector旋转到[-observe_dist, 0 , 0]
                    polar_vector_scan = polar_vector/rotate_factor_polar/rotate_factor;                 % 不同观测方向上的极化矢量
                    
                    [TB1(n,1),count1(n,1),e_tg1(n,1)] = FindFacet(inc_vector,observe_position,Param,polar_vector_scan);           % 计算像素点对应目标面元的亮温值
                    if TB1(n,1)==0               % 若像素点不在目标上，即入射射线与目标不相交，直接应用拟合系数计算面元亮温
                        zenith_theta = acosd(inc_vector*zenith_vector/norm(inc_vector));                % 计算入射射线与地表法向量的夹角，即为天顶角
                        zenith_theta = 90 - zenith_theta;
                        index = ceil((zenith_theta-(-90)+0.00001)/0.1);              % 判断俯仰角所在区间
                        factor = ceil((zenith_theta-(-90)+0.00001)/0.1)-(zenith_theta-(-90)+0.00001)/0.1;                        % 加权因子
                        TB1(n,1) = radio_atmos1(index)*factor+radio_atmos1(index+1)*(1-factor);        % 加权计算对应俯仰角下的背景辐射
                    end
                end
            end
           parfor_progress; 
        end
        TB(:,m) = TB1;
        e_tg(:,m) = e_tg1;
        count(:,m) = count1;
    end
    parfor_progress(0);
    
    parfor_progress(length(z_list));
    parfor m = 1:length(z_list)
        radio_atmos1 = radio_atmos;
        TB1 = TB(m,:);
        for n = 1:length(y_list)
            y_loc = observe_dist*tand(y_list(n));                   % 像素点在成像平面上的y坐标
            z_loc = tand(z_list(m))*observe_dist/cosd(y_list(n));      % 像素点在成像平面上的z坐标
            if TB1(1,n)==0       %若像素点不在目标区域内
                scan_vector = [-observe_dist, y_loc , z_loc];        % 旋转后的坐标系中入射射线的方向矢量
                inc_vector = scan_vector/rotate_factor;           % 入射射线在原始坐标系的方向矢量
                zenith_theta = acosd(inc_vector*zenith_vector/norm(inc_vector));                % 计算入射射线与地表法向量的夹角，即为天顶角
                zenith_theta = 90 - zenith_theta;
                index = ceil((zenith_theta-(-90)+0.00001)/0.1);              % 判断俯仰角所在区间
                factor = ceil((zenith_theta-(-90)+0.00001)/0.1)-(zenith_theta-(-90)+0.00001)/0.1;                        % 加权因子
                TB1(1,n) = radio_atmos1(index)*factor+radio_atmos1(index+1)*(1-factor);        % 加权计算对应俯仰角下的背景辐射
            end
        end
        TB(m,:) = TB1;
        parfor_progress;
    end
    parfor_progress(0);
    pixel_area = (observe_dist*tand(Param.Observe.pixel_space))^2;
    target_area = pixel_area*nnz(count);
    TB_target = sum(TB(count~=0))/nnz(count);                        % 目标区域平均亮温
end
%% 背景和目标一起计算
if Param.calculation_mode == 3                  % 并行计算
    observe_dist = Param.Observe.observe_dist*1000;  
%     polar_angle = Param.Observe.polar_angle;    % 极化角
    
     % 计算包围盒范围
    x_node = [max(trans_point(:,1)),min(trans_point(:,1))]-observe_dist;
    y_node = -[max(trans_point(:,2)),min(trans_point(:,2))]*observe_dist;
    z_node = -[max(trans_point(:,3)),min(trans_point(:,3))]*observe_dist;
    box_y = y_node'*(1./x_node);
    box_z = z_node'*(1./x_node);
    z_max = max(max(box_z));
    y_max = max(max(box_y));
    z_min = min(min(box_z));
    y_min = min(min(box_y));
    
    hwait = waitbar(0,'亮温计算');
%   for m=1:length(y_list)
%       for n = 1:length(z_list)
    for m = [38 78]           % 列
        for n = 87       % 行
            y_loc = observe_dist*tand(y_list(m));                   % 像素点在成像平面上的y坐标
            z_loc = tand(z_list(n))*observe_dist/cosd(y_list(m));      % 像素点在成像平面上的z坐标
            scan_vector = [-observe_dist, y_loc , z_loc];        % 旋转后的坐标系中扫描射线的方向矢量
            inc_vector = scan_vector/Param.Observe.rotate_factor;       % 扫描射线在原始坐标系的方向矢量
%             ph = cross(inc_vector,zenith_vector)/norm(cross(inc_vector,zenith_vector));
%             pv = cross(ph,inc_vector)/norm(cross(ph,inc_vector));
%             polar_vector_scan = ph*sind(polar_angle)+pv*cosd(polar_angle);         
            a = sind(z_list(n));
            b = cosd(z_list(n));
            c = sind(y_list(m));
            d = cosd(y_list(m));
            rotate_factor_polar = [b*d c -a*d
                -b*c d a*c
                a  0  b];                    % 坐标转换矩阵，由scan_vector旋转到[-observe_dist, 0 , 0]
            polar_vector_scan = Param.Observe.polar_vector/rotate_factor_polar/Param.Observe.rotate_factor;             % 不同观测方向上的极化矢量
            if y_loc<=y_max&&y_loc>=y_min&&...
                    z_loc<=z_max&&z_loc>=z_min
                [TB(n,m),count(n,m),e_tg(n,m)] = FindFacet(inc_vector,Param.Observe.observe_position,Param,polar_vector_scan);       % 计算像素点对应目标面元的亮温值
            end
            if TB(n,m)==0       %若像素点不在目标区域内或不在目标上
                zenith_theta = acosd(inc_vector*zenith_vector/norm(inc_vector));            % 反余弦求夹角（角度）
                zenith_theta = 90 - zenith_theta;
                index = ceil((zenith_theta-(-90)+0.00001)/0.1);              % 判断俯仰角所在区间
                factor = ceil((zenith_theta-(-90)+0.00001)/0.1)-(zenith_theta-(-90)+0.00001)/0.1;                        % 加权因子
                TB(n,m) = radio_atmos(index)*factor+radio_atmos(index+1)*(1-factor);        % 加权计算对应俯仰角下的背景辐射
            end
        end
        waitbar(m/length(y_list),hwait,'亮温计算');
    end
    close(hwait)
    pixel_area = (observe_dist*tand(Param.Observe.pixel_space))^2;
    target_area = pixel_area*nnz(count);
    TB_target = sum(TB(count~=0))/nnz(count);                        % 目标区域平均亮温
end
end
