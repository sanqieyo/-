function [TB,trace_count,e_tg] = RadiationTrace(facet_num,inc_vector,inc_node,trace_count,polar_vector,Param)

% 本函数用于复杂场景的亮温追踪
% 输入：面元编号facet_num，入射向量inc_vector，入射点inc_node，射线追踪次数trace_count，极化向量polar_vector
% 输出：反射方向亮温值TB，射线追踪次数trace_count
% rank_point = screening_rank_point;                   % 导入按面元顶点排列（x,y,z）坐标
% normal_vector = screening_normal_vector;             % 导入面元法向量
% target_rank = screening_target_rank;                 % 导入面元节点信息
rank_point = Param.Target.boat.rank_point;                   % 导入按面元顶点排列（x,y,z）坐标
normal_vector = Param.Target.boat.normal_vector;             % 导入面元法向量
target_rank = Param.Target.boat.target_rank;                 % 导入面元节点信息
target_point = Param.Target.boat.target_point;           % 导入节点坐标信息

[inc_theta,zenith_theta,ref_vector] = GetTheta(inc_vector,normal_vector(facet_num,:));              % 计算入射角、反射向量天顶角以及反射向量
[eh_tg,ev_tg] = RefCoherer(inc_theta,Param);         % 入射角为inc_theta，计算入射到飞机目标表面时发射率的水平极化和垂直极化分量
e_tg = LinearPolarization(polar_vector,inc_vector,normal_vector(facet_num,:),eh_tg,ev_tg);          % 由极化分量计算目标表面发射率
polar_vector = LinearPolarRota(polar_vector,inc_vector,ref_vector,normal_vector(facet_num,:));      % 旋转后的极化向量
T = Param.Target.boat.target_temp(facet_num);            % 目标面元温度


trace_count = trace_count+1;        % 记录射线追踪次数
if trace_count>50
    TB = e_tg*T+(1-e_tg)*GetEnvironRadiation(polar_vector,90-zenith_theta,Param,ref_vector);
    return;
end

facet_list1 = normal_vector*(ref_vector)';          % 通过入射向量与面元法向量的关系，判断该面元是否可见
facet_list2 = find(facet_list1<0);                  % 剔除不可见面元（若夹角小于targetdown_angle，则认为不可见）

%% 平面1求交运算

plane1_nvector = cross(ref_vector,inc_node);        % 平面1的法向量定义为入射点到坐标原点方向向量与入射向量的叉积
if plane1_nvector == zeros(1,3)         % 若平面1法向量为（0,0,0），将平面1法向量置为（1,0,0）
   plane1_nvector = [0 1 0];
end
plane1_nvector = plane1_nvector/max(abs(plane1_nvector));           % 简化平面1法向量的坐标，提高运算效率
plane1_result = plane1_nvector*target_point';           % 由点到面的距离公式，计算面元节点到平面1距离（只需判断正负，忽略分母）

Tr_a = sign(plane1_result(target_rank(facet_list2,1)));             % 面元第一个顶点到平面1距离的正负（正为1，负为-1，若在平面上则为0）
Tr_b = sign(plane1_result(target_rank(facet_list2,2)));             % 面元第二个顶点到平面1距离的正负（正为1，负为-1，若在平面上则为0）
Tr_c = sign(plane1_result(target_rank(facet_list2,3)));             % 面元第三个顶点到平面1距离的正负（正为1，负为-1，若在平面上则为0）
Tr_result = Tr_a+Tr_b+Tr_c;             % 判断面元顶点是否在平面1同一侧             
facet_list3 = facet_list2(Tr_result~=3&Tr_result~=-3);              % 选取面元顶点不在平面1同一侧的（全为正或全为负）
if isempty(facet_list3)         % 如果不存在相交面元，由入射向量直接计算面元亮温
    TB = e_tg*T+(1-e_tg)*GetEnvironRadiation(polar_vector,90-zenith_theta,Param,ref_vector);
    return
end

%% 平面2求交运算

plane2_nvector = cross(plane1_nvector,ref_vector);              % 平面2的法向量定义为平面1法向量与入射向量的叉积
if plane2_nvector == zeros(1,3)         % 若平面1法向量为（0,0,0），将平面1法向量置为（1,0,0）
   plane2_nvector = [0 0 1];
end
plane2_nvector = plane2_nvector/max(abs(plane2_nvector));       % 简化平面2法向量的坐标，提高运算效率
plane2_result=plane2_nvector*target_point';         % 由点到面的距离公式，计算面元节点到平面2距离Ax+By+Cz部分（只需判断正负，忽略分母）
d = plane2_nvector*(inc_node)';             % 计算D（-d）的值
plane2_result = plane2_result-d;            % 计算Ax+By+Cz+D的值，即计算面元节点到平面2距离

Tr_a = sign(plane2_result(target_rank(facet_list3,1)));             % 面元第一个顶点到平面2距离的正负（正为1，负为-1，若在平面上则为0）
Tr_b = sign(plane2_result(target_rank(facet_list3,2)));             % 面元第二个顶点到平面2距离的正负（正为1，负为-1，若在平面上则为0）
Tr_c = sign(plane2_result(target_rank(facet_list3,3)));             % 面元第三个顶点到平面2距离的正负（正为1，负为-1，若在平面上则为0）
Tr_result = Tr_a+Tr_b+Tr_c;             % 判断面元顶点是否在平面2同一侧
facet_list4 = facet_list3(Tr_result~=3&Tr_result~=-3);              % 选取面元顶点不在平面2同一侧的（全为正或全为负）
if isempty(facet_list4)         % 如果不存在相交面元，由入射向量直接计算面元亮温
    TB = e_tg*T+(1-e_tg)*GetEnvironRadiation(polar_vector,90-zenith_theta,Param,ref_vector);
    return
end

 %% 线面求交
 
 facet_list5 = zeros(size(facet_list4,1),1);        % 初始化与射线相交的候选面元列表
 re_node = zeros(size(facet_list4,1),3);            % 初始化与射线相交的候选面元交点
 counter = 0;           % 初始化与射线相交的候选面元个数
 
 for m = 1:size(facet_list4,1)  
     E_1 = rank_point(3*facet_list4(m)-1,:)-rank_point(3*facet_list4(m)-2,:);           % 三角面元顶点向量1
     E_2 = rank_point(3*facet_list4(m),:)-rank_point(3*facet_list4(m)-2,:);             % 三角面元顶点向量2
     V_0 = inc_node-rank_point(3*facet_list4(m,1)-2,:);             % 入射点到面元顶点向量
     D_0 = ref_vector/norm(ref_vector);             % 向量归一化
     node_result = [cross(V_0,E_1)*(E_2)';          % 由快速面元求交方法进行计算
             cross(D_0,E_2)*(V_0)';
             cross(V_0,E_1)*(D_0)']/(cross(D_0,E_2)*(E_1)');
    test2 = sign(node_result(1,1))+sign(node_result(2,1))+sign(node_result(3,1));       % 对计算结果变量进行判断  
    test1 = (node_result(2,1)+node_result(3,1))<1;
    if test2+test1==4           % 判断求解结果是否满足条件
        counter = counter+1;
       Re_Node = [1-(node_result(2,1)+node_result(3,1)),node_result(2,1),node_result(3,1)]* ...
            [rank_point(3*facet_list4(m)-2,:);rank_point(3*facet_list4(m)-1,:);rank_point(3*facet_list4(m),:)];   % 由计算结果解得射线与相交面元的交点坐标
       re_node(counter,:) = Re_Node(1,:);           % 得到射线与面元的交点
       facet_list5(counter,1) = facet_list4(m,1);   % 取得与射线相交的面元
    end
 end
 
 re_node = re_node(1:counter,:);                % 剔除不相交的面元
 facet_list5 = facet_list5(1:counter,:);        % 得到与平面2相交的面元列表
 if counter>0       % 判断交点是否存在  
     re_vector = re_node-repmat(inc_node,[counter,1]);      % 交点到入射点的方向向量
     node_dist = sum(abs(re_vector).^2,2).^(1/2);           % 计算各交点到入射点的距离
     facet_num = find(node_dist==min(node_dist));           % 找到距离入射点最近的交点
     facet_num = facet_num(1);                  % 避免交点出现在面元边上的情况
     [TB,trace_count,~] = RadiationTrace(facet_list5(facet_num),re_vector(facet_num,:),re_node(facet_num,:),trace_count,polar_vector,Param);
     TB = e_tg*T+(1-e_tg)*TB;
     % 对相交面元进行亮温追踪，返回亮温值
 else               % 不存在交点
     TB = e_tg*T+(1-e_tg)*GetEnvironRadiation(polar_vector,90-zenith_theta,Param,ref_vector);
 end
 
end