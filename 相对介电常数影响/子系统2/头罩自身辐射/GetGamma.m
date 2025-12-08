function [ Gamma_h, Gamma_v, L2, theta_t ] = GetGamma( Param, theta_in)

%% 计算天线罩表面反射率        
lambda = Param.SystemInput.lambda;  
d_flat = Param.radome.thickness;  % 天线罩厚度
lay_num = Param.radome.lay_num;       % 介质分层数
delta_d = d_flat/lay_num;

epsilon = Param.radome.epsilon; 
tan_delta = Param.radome.tan_delta; 
epsilon_c = epsilon * (1 - 1j * tan_delta);

alpha_2 = 2*pi/lambda*abs(imag(sqrt(epsilon_c)));
beta_2 = 2*pi/lambda*real(sqrt(epsilon_c));
k_alpha_2 = 2*alpha_2;  % 吸收系数

% 计算实透射角
k1 = 2*pi/lambda*sqrt(epsilon);
pp = 2*alpha_2*beta_2;
qq = beta_2^2 - alpha_2^2 - k1^2*(sin(theta_in)).^2;
theta_t = atan(sqrt(2)*k1*sin(theta_in)./(sqrt(sqrt(pp^2 + qq.^2) + qq)));

L2 = exp(k_alpha_2*delta_d./cos(theta_t));  % 损耗因子

% 计算功率反射系数，水平极化和垂直极化
Gamma_h = abs((cos(theta_in) - sqrt(epsilon_c - sin(theta_in).^2)) ./ (cos(theta_in) + sqrt(epsilon_c - sin(theta_in).^2))).^2;
Gamma_v = abs((epsilon_c * cos(theta_in) - sqrt(epsilon_c - sin(theta_in).^2)) ./ (epsilon_c * cos(theta_in) + sqrt(epsilon_c - sin(theta_in).^2))).^2;



end

