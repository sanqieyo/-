% 子系统2 不包含头罩幅相误差的理想可见度计算代码 2022.04.15
function [visibility_Ideal] = GetIdealVisibility( Param ,TB_original)

%% 输入场景参数处理
[ TB_scene, del_Fov_sc, ksai_sc, eta_sc, TB_origin, del_Fov_origin, ksai_origin, eta_origin] = ProcessTB( Param, TB_original);

%% 下一部分
UV_distribution = Param.Array.UV_distribution;
BandWidth = Param.SystemInput.BandWidth;
f0 = Param.SystemInput.f0;
Omega = Param.SystemInput.Omega;


%% 不考虑头罩影响的可见度计算（无噪声）
parfor mm = 1 : length(UV_distribution)  
    
       u = real(UV_distribution(mm));
       v = imag(UV_distribution(mm));
  
       fourier_vector = del_Fov_origin.*exp(-1j*2*pi*(u*ksai_origin+v*eta_origin)).*sinc(-BandWidth*(u*ksai_origin+v*eta_origin)/f0);
       visibility(mm) = TB_origin'*fourier_vector/Omega;             % 没有头罩影响的小视场可见度
             
       fourier_vector_large = del_Fov_sc.*exp(-1j*2*pi*(u*ksai_sc+v*eta_sc)).*sinc(-BandWidth*(u*ksai_sc+v*eta_sc)/f0);
       visibility_large(mm) = TB_scene'*fourier_vector_large/Omega;  % 没有头罩影响的背景可见度
       
end

visibility = visibility + visibility_large;  %没有头罩影响的可见度

visibility_Ideal = visibility.';

end

