function target_temp = ProcessTemp(target_speed)
% 根据目标速度计算目标面元温度
% 输入：目标飞行速度target_speed
% 输出：目标面元温度target_temp

column = ceil((target_speed-0.5+0.00001)/0.1);              % 判断飞行速度所在区间
factor = ceil((target_speed-0.5+0.00001)/0.1)-(target_speed-0.5+0.00001)/0.1;                        % 加权因子
target_temp = load('模型数据\f22_temp.txt');        % 导入目标面元温度信息
target_temp = target_temp(:,column)*factor+target_temp(:,column+1)*(1-factor);        % 加权计算对应速度下的目标面元温度    
end