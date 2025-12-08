function polar_vector=LinearPolarRota(polar_vector,inc_vector,ref_vector,normal_vector)
%函数功能：计算入射到平面时，极化角的变换
%输入：入射前极化向量polar_vector，入射向量inct，反射向量r_vector，面法向量Nt
%输出：反射后极化向量polar_vector

%入射前
ph1=cross(inc_vector,normal_vector)/norm(cross(inc_vector,normal_vector));%水平极化方向单位矢量
pv1=cross(ph1,inc_vector)/norm(cross(ph1,inc_vector));%垂直极化方向单位矢量

%反射后
ph2=cross(ref_vector,normal_vector)/norm(cross(ref_vector,normal_vector));%水平极化方向单位矢量
pv2=cross(ph2,ref_vector)/norm(cross(ph2,ref_vector));%垂直极化方向单位矢量

a = pv1*polar_vector'/(norm(pv1)*norm(polar_vector));     % 极化方向向量与pv的夹角的余弦
alpha = acosd(a);
alpha = real(alpha);

polar_vector = ph2*sind(alpha)+pv2*cosd(alpha);

end

