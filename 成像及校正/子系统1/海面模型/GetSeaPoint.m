function sea_point = GetSeaPoint(Param,tt)


%----------------海面高度生成--------------------%
height = oceanprodct(Param,tt);

%----------------海面高度矩阵生成--------------------%
n = length(height(1,:));
m = length(height(:,1));
sea_point = zeros(n * m, 3);
count = 1;
for ii = 1 : n
    for jj = 1 : m
        sea_point(count, 1) = ii * Param.delta_x;
        sea_point(count, 2) = jj * Param.delta_y;
        sea_point(count, 3) = height(ii, jj);
        count = count + 1;
    end
end        