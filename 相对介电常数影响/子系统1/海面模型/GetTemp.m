function temp_data = GetTemp(Param)

n = Param.x_length / Param.delta_x;
m = Param.y_length / Param.delta_y;
count = 2 * (n - 1) * (m - 1);


temp_data = zeros(count,1);
temp_data(:,1)=293;
