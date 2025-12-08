function [ TB_patch ] = GetPolarTB( Param, Gamma, L2, T_out, T_in )
% 根据不同极化反射率，计算相应的亮温
lay_num = Param.radome.lay_num;
N = lay_num+1;

T_layer = GetProfileT( T_out, T_in, lay_num );
Ts_layer = [0 (1-1/L2)*T_layer];

% 计算每一层的辐射亮温
Gamma_layer = zeros(1,lay_num+1); % 反射率
Gamma_layer(1) = Gamma;
Gamma_layer(length(Gamma_layer)) = Gamma;
L = L2*ones(1,lay_num+1);
L(1) = 0;
TB_temp = zeros(1,N);
for ii = 2 : N
    temp1 = 1;
    for j = 2 : ii
        temp1 = temp1 * (1-Gamma_layer(j-1))/(1-Gamma_layer(j-1)*Gamma_layer(j)/L(j)/L(j));
    end

    temp2 = 1;
    for k = 2 : ii-1
        temp2 = temp2 / L(k);
    end

    TB_temp(ii) = (1 + Gamma_layer(ii)/L(ii)) * Ts_layer(ii) * temp1 * temp2;
end

TB_patch = sum(TB_temp);

end

