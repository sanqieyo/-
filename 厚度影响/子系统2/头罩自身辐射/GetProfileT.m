% 此函数根据给定头罩内外表面温度和分层数，给出剖面温度分布（稳态均匀分布模型）
% T_profile表示剖面温度
function [ T_profile ] = GetProfileT( T_out, T_in, layer_num )

delta_T = (T_out-T_in)/(layer_num-1);

T_profile = zeros(1,layer_num);
for i = 1:layer_num
    T_profile(i) = T_in+(i-1)*delta_T;
end

end
