function [inc_theta,zenith_theta,ref_vector] = GetTheta(inc_vector,normal_vector)
% 计算当前面与入射方向的夹角，反射方向与法向量的夹角
% 输入：当前面法向量normal_vector、入射向量inc_vector  (行向量)
% 输出：当前面与探测方向的夹角（入射角）inc_theta,探测方向的反射方向与铅垂方向的夹角y2

inc_vector = inc_vector/norm(inc_vector);
normal_vector = normal_vector/norm(normal_vector);    % 将入射向量与法向量单位化
x1 = inc_vector*normal_vector';
inc_theta = 180-acosd(x1);         % 反余弦求夹角（角度）,原入射方向与面法向量夹角为钝角，此处得到的夹角与计算发射率时夹角一致
%y1
vector_v = inc_vector-(normal_vector*x1);   % 分解出与当前面面向量相垂直的向量vector_v,  (n_v.*x1)为入射向量沿法线方向的分量
ref_vector = -x1*normal_vector+vector_v;      % 得到入射向量经当前面反射后的方向向量
zenith_vector = [0;0;1];        % 铅垂方向向量
x2 = ref_vector*zenith_vector/norm(ref_vector);     % 夹角的余弦值
zenith_theta = acosd(x2);     % 反余弦求夹角（角度）

end