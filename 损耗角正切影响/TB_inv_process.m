function[TB] = TB_inv_process(ksai_INV,eta_INV,ksai_inv,eta_inv,in_blur,UV_distribution,del_uv,Omega,visibility)
% ksai_INV = linspace(-sind(5), sind(5), 1001);
% eta_INV = linspace(-sind(5), sind(5), 1001);
% [ksai_inv, eta_inv] = meshgrid(ksai_INV, eta_INV);
TB = zeros(length(ksai_INV),length(eta_INV)); % 无冗余基线反演

for mm = 1 : size(ksai_inv,1)
    for nn = 1 : size(ksai_inv,2)
        Ksai = ksai_inv(mm,nn);
        Eta = eta_inv(mm,nn);
        ifourier_vector = exp(1j*2*pi*(real(UV_distribution)*Ksai+imag(UV_distribution)*Eta));
        TB(mm,nn) = del_uv*Omega*visibility.'*ifourier_vector;
        %         a = visibility_FF.';
    end
end
TB = real(TB);
TB(~in_blur) = nan;


% figure;
% imagesc(ksai_INV, eta_INV, TB);
% colormap(jet);
% colorbar;
% shading interp;
% title('头罩影响的可见度--3');


% 让六边形视场外的TB值不显示
% TB = TB.*in_blur;
% [a, b] = find(TB==0);
% for i = 1:length(a)
%     TB(a(i),b(i)) = nan;
% end
% TB;
end