function [eh,ev]=RefGround(inc_theta,Param)
%计算地面的发射率
%输入：入射方向与垂直方向夹角，即入射角inc_theta
%输出：输出发射率的水平极化eh和垂直极化分量ev

jiedian = Param.ground_permittivity;
inc_theta = inc_theta*pi/180;%计算垂直入射角
Rh = abs((cos(inc_theta)-sqrt(jiedian-(sin(inc_theta)).^2))./(cos(inc_theta)+sqrt(jiedian-(sin(inc_theta)).^2))).^2;
Rv = abs((sqrt(jiedian-(sin(inc_theta)).^2)-jiedian.*cos(inc_theta))./(sqrt(jiedian-(sin(inc_theta)).^2)+jiedian.*cos(inc_theta))).^2;
eh = 1-Rh;
ev = 1-Rv;
end
