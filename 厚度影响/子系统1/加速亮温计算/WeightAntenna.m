function [antenna_temp,counter] = WeightAntenna( z_loc,y_loc,Param,position )
% 本函数用于计算各像素点处的天线温度
% 输入天线温度成像平面像素点坐标（z_loc,y_loc），亮温分布TB
% 输出天线温度antenna_temp，波束范围内亮温点个数counter

D = Param.observe_dist*1000;            % 观测距离
beita = Param.beam_width*pi/180;       % 天线波束宽度

d_wave = 0.3/Param.center_freq;         % 波长

counter = 0;
Dim = d_wave/beita;         % 天线口径
z_d = D/cosd(y_loc)*tand(z_loc);        % 扫描射线在投影平面上的纵坐标
y_d = D*tand(y_loc);        % 扫描射线在投影平面上的横坐标
d = sqrt(z_d^2+y_d^2);
r_max = tan(atan(d/D)+beita/2)*D-d;     % 圆形波束在投影平面上的最大半径

antenna_temp = 0;           % 天线温度初始化
scan_vector = [-D,y_d,z_d];         % 旋转后坐标系中天线波束中心点扫描射线
scan_vector = scan_vector/norm(scan_vector);        % 将扫描射线归一化
%% 矩阵计算法

% position = Param.pixel_position;        % 亮温分布各点坐标
position = position(abs(position(:,1)-y_d)<=r_max&abs(position(:,2)-z_d)<=r_max,:);         % 选取在波束最大半径范围内的点
inc_vector = [(-D)*ones(size(position,1),1),position(:,1),position(:,2)];               % 计算各点到观测点的入射向量
inc_vector = inc_vector./repmat(sum(abs(inc_vector).^2,2).^(1/2),[1,3]);                % 将入射向量inc_vector归一化
theta = acos(scan_vector*inc_vector')';         % 计算入射向量与扫描射线的夹角
logical = abs(theta)<=beita/2;          % 得到夹角小于1/2波束宽度的逻辑索引
theta = theta(logical,1);               % 由逻辑索引得到满足条件的theta值
TB = position(logical,3);               % 由逻辑索引得到满足条件的TB值
if ~isempty(theta)
    theta(theta==0)=eps;                % 若theta为0，则将其置为无穷小（eps）
    Fn = (2*besselj(1,pi*Dim*sin(theta)/d_wave)./(pi*Dim*sin(theta)/d_wave)).^2;          % 计算归一化理想功率方向图权值
    antenna_temp = sum(TB.*Fn);            % 天线功率方向图权值与亮温值加权
    F = sum(Fn);
    antenna_temp = antenna_temp/F;
    antenna_temp = real(antenna_temp);          % 由于由反余弦acos得到theta值过程中，存在误差，导致antenna_temp出现虚部，此处只取实部
    counter = length(theta);
end

end

% %% for循环计算法
% step = D*tand(Param.scan_step)*Param.scan_factor;        % 扫描角度步长
% y = Param.horizon_scale;
% z = Param.vertical_scale;
% d_3 = zeros(ceil(r_max/step)^2,3);          % 用于存放可能位于波束范围内的亮温点坐标
% d_5 = zeros(1,ceil(r_max/step)^2);          % 用于存放可能位于波束范围内的亮温点的值
% r_count = 0;             % 初始化波束范围内的亮温点个数
% for m = 1:size(z,2)      % 取得可能在波束范围内的亮温点
%     if abs(z(m)-z_d)<=r_max
%         for n=1:size(y,2)
%             if abs(y(n)-y_d)<=r_max
%                 r_count=r_count+1;
%                 d_2=[-D,y(n),z(m)];       % 旋转后坐标系中亮温点与观测点间的方向矢量
%                 d_2=d_2/norm(d_2);        % 将方向矢量归一化
%                 d_3(r_count,:)=d_2;
%                 d_5(r_count)=TB(m,n);
%             end
%         end
%     end
% end
% d_3 = d_3(1:r_count,:);         % 截取包含亮温点部分
% d_5 = d_5(:,1:r_count);
% theta1 = acos(abs(scan_vector*d_3'));         % 亮温点矢量与波束中心点的夹角
% d_6 = find(theta1<=beita&theta1>0);           % 判断夹角是否小于波束宽度
% d_7 = theta1(d_6);               % 截取波束范围内点与中心点d的夹角
% d_8 = d_5(d_6);                 % 截取波束范围内点的亮温值
% 
% if ~isempty(d_6)
%     Fn = (2*besselj(1,pi*Dim*sin(d_7)/d_wave)./(pi*Dim*sin(d_7)/d_wave)).^2;              % 计算归一化理想功率方向图权值
%     antenna_temp = sum(d_8.*Fn);              % 天线功率方向图权值与亮温值加权
%     F = sum(Fn);
%     antenna_temp = antenna_temp/F;
%     counter = length(d_6);
% end

% end


