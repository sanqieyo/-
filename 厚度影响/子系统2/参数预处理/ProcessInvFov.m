function [ ksai_INV, eta_INV, ksai_inv, eta_inv, ksai_inv_rec, eta_inv_rec, in_blur ] = ProcessInvFov( Param, del_v )
%UNTITLED9 此处显示有关此函数的摘要
%   此处显示详细说明
inv_sc_pix_count = Param.inv_sc_pix_count;
Scene_Xi_inv = Param.Scene_Xi_inv;
Scene_Eta_inv = Param.Scene_Eta_inv;
ksai_INV=linspace(-Scene_Xi_inv,Scene_Xi_inv,inv_sc_pix_count);%!!!!!!场景的像素点大小必须小于系统主波束的1/3
eta_INV=linspace(-Scene_Eta_inv,Scene_Eta_inv,inv_sc_pix_count);%此处修改像素数
[ksai_inv,eta_inv]=ndgrid(ksai_INV,eta_INV);
% [ksai_inv,eta_inv]=meshgrid(ksai_INV,eta_INV);
ksai_inv_rec = reshape(ksai_inv,[],1);
eta_inv_rec = reshape(eta_inv,[],1);

del_u_new = del_v;
x = linspace(0, 2*pi, 7);
X = exp(1j*x);
a = 1j*linspace(-1/del_u_new/2/sqrt(3),(1/del_u_new/2/sqrt(3)),20)+1/del_u_new/2;
fov2 = [a a*X(1) a*X(2) a*X(3) a*X(4) a*X(5) a*X(6)];
fov = fov2;
fov = fov*exp(1j*2*pi/12);
in_plot_figure = inpolygon(ksai_inv_rec,eta_inv_rec,real(fov),imag(fov));
in_blur = reshape(in_plot_figure,length(ksai_INV),length(eta_INV));

% figure
% plot(real(fov),imag(fov),eta_inv_rec(in_blur),ksai_inv_rec(in_blur));
end
