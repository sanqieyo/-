function [ rank_point, normal_vector ] = ProcessRadome( Param )

%UNTITLED10 此处显示有关此函数的摘要
%   此处显示详细说明
% 输出面元各顶点的坐标rank_point，面元法向量normal_vector
patch_matrix = Param.patch_matrix;
point_matrix = Param.point_matrix;
% patch_matrix = Param.patch_matrix_part;
% point_matrix = Param.point_matrix_part;

rank_point = zeros(3*size(patch_matrix.patch_matrix,1),3);   %存储按面元排列的结点坐标
for m=1:size(patch_matrix.patch_matrix,1)
    rank_point(3*m-2,:)=point_matrix.point_matrix(patch_matrix.patch_matrix(m,1),1:3);
    rank_point(3*m-1,:)=point_matrix.point_matrix(patch_matrix.patch_matrix(m,2),1:3);
    rank_point(3*m,:)=point_matrix.point_matrix(patch_matrix.patch_matrix(m,3),1:3);
end

normal_vector=zeros(size(patch_matrix.patch_matrix,1),3);        %法向量矩阵初始化
 %incident_point=zeros(size(Param.target_rank,1),3);       %入射点矩阵初始化
for n=1:size(patch_matrix.patch_matrix,1)
    a1=rank_point(3*n-1,:)-rank_point(3*n-2,:);
    a2=rank_point(3*n,:)-rank_point(3*n-1,:);       %计算三角面元各边的方向矢量
    normal_vector(n,:)=cross(a1,a2);                %三角面元两条边的叉积即为法向量
         %incident_point(n,:)=(rank_point(3*n,:)+rank_point(3*n-1,:)+rank_point(3*n-2,:))/3;      %面元三个顶点的中心点
end


end

