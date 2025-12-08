function polar_vector=ProcessPolarization(Param)
%本函数用于计算极化向量
%输入极化角polar_angle，坐标系旋转矩阵rotate_factor
%输出极化向量polar_vector
ph = [0,1,0];
pv = [0,0,1];
polar_vector = ph*sind(Param.polar_angle)+pv*cosd(Param.polar_angle);
end