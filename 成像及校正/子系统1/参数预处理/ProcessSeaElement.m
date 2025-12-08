function [rank_point,normal_vector,target_point,shift_x,shift_y,shift_z]=ProcessSeaElement(Param,Observe, a, b)
%本函数用于计算由面元各顶点坐标组成的坐标矩阵
%输入结构体SimuParam，目标高度target_altitude
%输出面元各顶点的坐标rank_point，面元法向量normal_vector，面元入射点（即中心点）incident_point
target_point = Param.sea.target_point;
target_rank = Param.sea.target_rank;
target_altitude = Observe.target_altitude;
observe_position = Observe.observe_position_1;

shift_x = (max(target_point(:,1))+min(target_point(:,1)))/2;
target_point(:,1) = target_point(:,1) - shift_x;      
% 将节点x坐标中心移至目标中心
shift_y = (max(target_point(:,2))+min(target_point(:,2)))/2; 
target_point(:,2) = target_point(:,2) - shift_y; 
% 将节点y坐标中心移至目标中心
shift_z = 0;
if target_altitude>(max(target_point(:,3))-min(target_point(:,3)))/2
    shift_z = (max(target_point(:,3)) + min(target_point(:,3)))/2;
    target_point(:,3) = target_point(:,3) - shift_z;
end
% 将节点z坐标中心移至目标中心 
roll_angle = Param.roll_angle;          % x轴旋转角度，即翻滚角
pitch_angle = Param.pitch_angle;        % y轴旋转角度，即俯仰角
yaw_angle = Param.yaw_angle;            % z轴旋转角度，即偏航角
rotate_roll = [1 0 0;0 cosd(roll_angle) sind(roll_angle);0 -sind(roll_angle) cosd(roll_angle)];            % 翻滚角旋转矩阵
rotate_yaw = [cosd(yaw_angle) sind(yaw_angle) 0;-sind(yaw_angle) cosd(yaw_angle) 0;0 0 1];            % 偏航角旋转矩阵

rotate_axis = cross([1 0 0]*rotate_yaw,[0 0 1]);                    % 俯仰角旋转轴
rotate_axis = rotate_axis/norm(rotate_axis);                        % 转换为单位矢量
n1 = rotate_axis(1);                    % 旋转轴单位矢量x方向分量
n2 = rotate_axis(2);                    % 旋转轴单位矢量y方向分量
n3 = rotate_axis(3);                    % 旋转轴单位矢量z方向分量
rotate_pitch = [n1^2*(1-cosd(pitch_angle))+cosd(pitch_angle),n1*n2*(1-cosd(pitch_angle))+n3*sind(pitch_angle),...
    n1*n3*(1-cosd(pitch_angle))-n2*sind(pitch_angle);n1*n2*(1-cosd(pitch_angle))-n3*sind(pitch_angle),...
    n2^2*(1-cosd(pitch_angle))+cosd(pitch_angle),n2*n3*(1-cosd(pitch_angle))+n1*sind(pitch_angle);...
    n1*n3*(1-cosd(pitch_angle))+n2*sind(pitch_angle),n2*n3*(1-cosd(pitch_angle))-n1*sind(pitch_angle),...
    n3^2*(1-cosd(pitch_angle))+cosd(pitch_angle)];            % 俯仰角旋转矩阵
target_point = target_point*rotate_roll*rotate_yaw*rotate_pitch;                 % 将目标坐标旋转至不同飞行姿态

%%      考虑导弹头移动对成像过程的影响
theta = -atan(observe_position(1,2) / observe_position(1,1));
rotate_factor = [cos(theta) sin(theta) 0
                 -sin(theta) cos(theta) 0
                 0 0 1];
observe_position = observe_position * rotate_factor - b;
observe_position = observe_position / rotate_factor;

rotate_factor_1 = Observe.rotate_factor;
inc_vector = [a(1,1)*1000 0 0];
inc_vector = inc_vector / rotate_factor_1;
x_0 = observe_position(1);
y_0 = observe_position(2);
z_0 = observe_position(3);
x_loc = x_0 - inc_vector(1) / inc_vector(3) * z_0; 
y_loc = y_0 - inc_vector(2) / inc_vector(3) * z_0;

target_point(:,1) = target_point(:,1) - x_loc;
target_point(:,2) = target_point(:,2) - y_loc;

shift_x = shift_x + x_loc;
shift_y = shift_y + y_loc;

%%
rank_point = zeros(3*size(target_rank,1),3);   %存储按面元排列的结点坐标
for m=1:size(target_rank,1)
    rank_point(3*m-2,:)=target_point(target_rank(m,1),1:3);
    rank_point(3*m-1,:)=target_point(target_rank(m,2),1:3);
    rank_point(3*m,:)=target_point(target_rank(m,3),1:3);
end

normal_vector=zeros(size(target_rank,1),3);        %法向量矩阵初始化
 %incident_point=zeros(size(Param.target_rank,1),3);       %入射点矩阵初始化
for n=1:size(target_rank,1)
    a1=rank_point(3*n-1,:)-rank_point(3*n-2,:);
    a2=rank_point(3*n,:)-rank_point(3*n-1,:);       %计算三角面元各边的方向矢量
    normal_vector(n,:)=cross(a1,a2);                %三角面元两条边的叉积即为法向量
         %incident_point(n,:)=(rank_point(3*n,:)+rank_point(3*n-1,:)+rank_point(3*n-2,:))/3;      %面元三个顶点的中心点
end

end

