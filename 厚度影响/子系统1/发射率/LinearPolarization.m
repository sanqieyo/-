function e_tg  =  LinearPolarization(polar_vector,inc_vector,normal_vector,eh_tg,ev_tg)
% 根据相应的极化角计算任意线极化角的辐射率
% 输入：极化方向向量polar_vector、入射向量inc_vector、面法向量normal_vector、水平极化辐射率eh_tg、垂直极化辐射率ev_tg
% 输出：对应面的线极化辐射率e_tg，极化向量polar_vector

% 入射前，求解h、v轴向量
if(cross(inc_vector,normal_vector) == 0)        % 垂直入射时
    e_tg = eh_tg;
else
    ph = cross(inc_vector,normal_vector)/norm(cross(inc_vector,normal_vector));     % 水平极化方向单位矢量
    pv = cross(ph,inc_vector)/norm(cross(inc_vector,ph));   % 垂直极化方向单位矢量
    
    %求解极化方向向量与pv的夹角
    a = pv*polar_vector'/(norm(pv)*norm(polar_vector));     % 极化方向向量与pv的夹角的余弦
    alpha = acosd(a);
    alpha = real(alpha);
    
    e_tg = eh_tg*cosd(alpha)^2+ev_tg*sind(alpha)^2;           % 极化组合到相应的极化角的线极化辐射率
end