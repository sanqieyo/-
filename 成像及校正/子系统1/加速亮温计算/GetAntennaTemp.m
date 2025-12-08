function antenna_temp = GetAntennaTemp( Param,TB )
% 本函数用于天线温度的加权计算
%   Detailed explanation goes here
% tic
horizon_scale = Param.horizon_scale;                % 亮温分布成像平面横向序列
vertical_scale = Param.vertical_scale;              % 亮温分布成像平面纵向序列
observe_dist = Param.observe_dist*1e3;              % 观测距离，m

pixel_position = zeros(length(horizon_scale)*length(vertical_scale),3);     % 将亮温平面像素点各点位置坐标初始化
for m = 1:length(horizon_scale)
    a = 1+(m-1)*length(vertical_scale);
    b = a+length(vertical_scale)-1;
    pixel_position(a:b,1) = ones(length(vertical_scale),1)*observe_dist*tand(horizon_scale(m));         % 第1列存储y坐标
    pixel_position(a:b,2) = observe_dist/cosd(horizon_scale(m))*tand(vertical_scale');                  % 第2列存储z坐标
%     pixel_position(a:b,4) = ones(length(vertical_scale),1)*horizon_scale(m);      % 第4列存储所在列角度 
%     pixel_position(a:b,5) = (vertical_scale)';         % 第5列存储所在行角度 
    pixel_position(a:b,3) = TB(:,m);                % 第3列存储亮温值
end

% beam_width = Param.beam_width/2;
Ymax = Param.scan_step*Param.horizon_scan/2;        % 由扫描次数确定视场范围
Zmax = Param.scan_step*Param.vertical_scan/2;      
y = linspace(-Ymax,Ymax,Param.horizon_scan);        % 建立天线温度成像平面
z = linspace(Zmax,-Zmax,Param.vertical_scan); 
antenna_temp = zeros(length(z),length(y));          % 天线温度初始化
% hwait = waitbar(0,'天线方向图');
% for m = 1:length(z)
%     for n = 1:length(y)
%         [antenna_temp(m,n),~] = WeightAntenna(z(m),y(n),TB,Param);
%     end
%     waitbar(m/length(z),hwait,'天线方向图');
% end
% close(hwait)
y_length = length(y);
parfor_progress(length(z));
parfor m = 1:length(z)
    TB1 = zeros(1,y_length);
   for n = 1:length(y)
       [TB1(1,n),~] = WeightAntenna(z(m),y(n),Param,pixel_position);
   end
    antenna_temp(m,:)=TB1;
    parfor_progress;
end
parfor_progress(0);
% toc
end

