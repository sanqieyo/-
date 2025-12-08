function sea_rank = GetSeaRank(Param)
%--------------按逆时针方向获得对应大小的rank矩阵-----------------%


n = Param.x_length / Param.delta_x;
m = Param.y_length / Param.delta_y;
count = 2 * (m - 1) * (n - 1);
sea_rank = ones(count, 3);
for ii = 0 : count - 1
    rowIndex = fix(ii / (2 * m - 2));
    cluIndex1 = fix(mod(ii, (2 * m - 2)) / 2);
    cluIndex2 = mod(mod(ii , 2 * m - 2), 2);
    sea_rank(ii + 1, 1) = rowIndex * m + cluIndex1 + 1;
    if cluIndex2 == 0
        sea_rank(ii + 1, 2) = (rowIndex + 1) * m + cluIndex1 + 1;
        sea_rank(ii + 1, 3) = (rowIndex + 1) * m + cluIndex1 + 2;
    else
        sea_rank(ii + 1, 2) = (rowIndex + 1) * m + cluIndex1 + 2;
        sea_rank(ii + 1, 3) = rowIndex * m + cluIndex1 + 2;
    end   
end 