function [ screening_target_rank, screening_normal_vector, screening_rank_point] = PreScreening( inc_vector, observe_position, Param )
%PreScreening 利用此函数预先筛选出第一个入射点十米内的面元

%------------------计算入射光线在原始坐标系中的具体位置---------%

x_0 = observe_position(1) + Param.Target.shift_x; 
y_0 = observe_position(2) + Param.Target.shift_y; 
z_0 = observe_position(3) + Param.Target.shift_z;

%------------------计算像素点射线在目标坐标系中的位置-----------%
x_loc = x_0 - inc_vector(1) / inc_vector(3) * z_0; 
y_loc = y_0 - inc_vector(2) / inc_vector(3) * z_0;

delta_x = Param.Target.delta_x;
delta_y = Param.Target.delta_y;

rowIndex = round(x_loc / delta_x);
cluIndex = round(y_loc / delta_y);

target_point = Param.Target.sea.target_point;
target_rank = Param.Target.sea.target_rank;
normal_vector = Param.Target.sea.normal_vector;  
n = round(Param.Target.x_length / delta_x) - 1;

rowIndex = max(1, rowIndex - 20);
cluIndex = max(1, cluIndex - 20);

rowIndex = min(n - 39, rowIndex);
cluIndex = min(n - 39, cluIndex);

%---------------------------更新线面求交的三个面元矩阵----------%
k = 1;
for ii = rowIndex + 1 : rowIndex + 39
    for jj = cluIndex + 1 : cluIndex + 39
        screening_target_rank(k,:) = target_rank((ii - 1) * 2 * n + 2 * jj - 1,:);
        screening_normal_vector(k,:) = normal_vector((ii - 1) * 2 * n + 2 * jj - 1,:);
        k = k + 1;
        screening_target_rank(k,:) = target_rank((ii - 1) * 2 * n + 2 * jj ,:);
        screening_normal_vector(k,:) = normal_vector((ii - 1) * 2 * n + 2 * jj,:);
        k = k + 1;
    end    
end
% k = 1;
% for ii = 1 : 99
%     for jj = 1 : 99
%         screening_target_rank(k,:) = (ii - 1) * 2 * (n - 1) + 2 * jj -   1;
%         screening_normal_vector(k,:) = (ii - 1) * 2 * (n - 1) + 2 * jj - 1;
%         k = k + 1;
%         screening_target_rank(k,:) = (ii - 1) * 2 * (n - 1) + 2 * jj;
%         screening_normal_vector(k,:) = (ii - 1) * 2 * (n - 1) + 2 * jj;
%         k = k + 1;
%     end    
% end

screening_rank_point = zeros(3*size(screening_target_rank,1),3);   %存储按面元排列的结点坐标
for m=1:size(screening_target_rank,1)
    screening_rank_point(3*m-2,:)=target_point(screening_target_rank(m,1),1:3);
    screening_rank_point(3*m-1,:)=target_point(screening_target_rank(m,2),1:3);
    screening_rank_point(3*m,:)=target_point(screening_target_rank(m,3),1:3);
end

end

