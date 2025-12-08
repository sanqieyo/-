function y=itu676d(f,P,T,ro)
%
% Group delay in picoseconds/km due to the dispersion of atmospheric gases
%
% f= frequency (GHz)
% P = total pressure = p+e in hpa
% T = temperature in Kelvin
% ro = water vapour concentration (g/m^3)
% p = dry air pressure in hPa
% e = water-vapour partial pressure in hpa
%
% DRAFT REVISION OF RECOMMENDATION ITU-R P.676-5 (29/11/2004)
%
% Created by Xiong Zubiao, 02/04/2005
%
%
global fw b1 b2 b3 b4 b5 b6 fo a1 a2 a3 a4 a5 a6
if isempty(fw)==1 | fw==0,
load wat676
load oxy676
end
%
theta=300.0./T;
e=ro.*T./216.7;
p=P-e;
%
% Line strength Si
%
Sio=a1.*(1e-7).*p.*(theta.^3.0).*exp(a2.*(1-theta));
Siw=b1.*(1e-1).*e.*(theta.^3.5).*exp(b2.*(1-theta));
%
% Line shape factor
%
deltao=(a5+a6.*theta).*(1e-4).*p.*theta.^(0.8);
deltaw=0;
deltafo=a3.*1e-4.*(p.*theta.^(0.8-a4)+1.1.*e.*theta);
deltafw=b3.*1e-4.*(p.*theta.^b4+b5.*e.*theta.^b6);
%
%
part1Fio=((fo-f)+deltafo.*(deltafo+f.*deltao)./fo)./((fo-f).^2+deltafo.^2);
part2Fio=((fo+f)+deltafo.*(deltafo-f.*deltao)./fo)./((fo+f).^2+deltafo.^2);
%
Fio=(part1Fio+part2Fio)-(2.0./fo);
%
%
%   Water vapour
%
part1Fiw=((fw-f)+deltafw.*(deltafw+f.*deltaw)./fw)./((fw-f).^2+deltafw.^2);
part2Fiw=((fw+f)+deltafw.*(deltafw-f.*deltaw)./fw)./((fw+f).^2+deltafw.^2);
%
Fiw=(part1Fiw+part2Fiw)-(2.0./fw);
%
%
% Dry air continuum N1Df ( N'D(f) )
%
d=4.8e-4.*(p+1.1.*e).*theta;
%
N1Df=6.14e-5.*p*theta.^2.0.*((1.0./(1+(f./d).^2))-1);
%
%
% N1Wf=0.998.*f.^2.0.*(1-0.20.*theta).*1.0e-6.*e.*theta.^2.7;
%
N1f=sum(Fio.*Sio)+sum(Fiw.*Siw)+N1Df;
%
y=N1f;
