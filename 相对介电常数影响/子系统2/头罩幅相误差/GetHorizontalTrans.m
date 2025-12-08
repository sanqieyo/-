% 该函数计算头罩的水平极化传输系数和相位延迟，输入theta单位为°
function [ trans_h, phase_delay_h ] = GetHorizontalTrans( theta,d ,epsilon,delta,f0)
%% 基本参数设置
com_epsilon = epsilon * (1 - 1j * delta);
theta_rad = theta * pi /180;
lambda = 3e8/f0;

com_admittance = com_epsilon*cos(theta_rad)/(sqrt(com_epsilon-(sin(theta_rad)^2))); %水平极化复导纳
com_phase = (2*pi*d/lambda)*sqrt(com_epsilon-(sin(theta_rad)^2));  % 复相角

%% 计算传输矩阵
trans_matrix(1,1) = cos(com_phase);
trans_matrix(1,2) = 1j*sin(com_phase)/com_admittance;
trans_matrix(2,1) = 1j*sin(com_phase)*com_admittance;
trans_matrix(2,2) = cos(com_phase);

E_H_incident = trans_matrix*[1;1];
T = 2/(E_H_incident(1,1)+E_H_incident(2,1)); 
trans_h = abs(T); % 电压传输系数

T_angle = angle(T)./pi*180;
for i = 1 : length(T_angle)
    if T_angle(i) < 0
        T_angle(i) = T_angle(i) + 360;
    end
end
delta_phi = 360/lambda*d.*cos(theta_rad);
while delta_phi > 360
    delta_phi = delta_phi - 360;
end
while delta_phi < 0
    delta_phi = delta_phi + 360;
end

% 平行于介质分界面上的传输系数相位分量
theta_1 = asin(sin(theta_rad)./sqrt(epsilon)); % 折射角
T_angle_horizontal = 2*pi/lambda*sin(theta_rad).*d.*tan(theta_1);
T_angle_horizontal = T_angle_horizontal./pi*180;
while T_angle_horizontal < 0
    T_angle_horizontal = T_angle_horizontal + 360;
end
while T_angle_horizontal > 360
    T_angle_horizontal = T_angle_horizontal - 360;
end

% IPD = -T_angle - delta_phi; 
IPD = -T_angle + T_angle_horizontal - delta_phi;  % 加上平行于介质分界面上的传输系数相位分量
while IPD < 0
    IPD = IPD + 360;
end
while IPD > 360
    IPD = IPD - 360;
end

phase_delay_h = IPD;
phase_delay_h = phase_delay_h*pi/180;
