function y=itu676a(f,P,T,ro)

% Absorption in dB/km due to atmospheric gases
%
% f = frequency (GHz)
% P = total pressure = p+e in hPa
% T = temperature in Kelvin
% ro = water vapour concentration (g/m^3)
% p = dry air pressure in hPa
% e = water-vapour partial pressure 
%
% DRAFT REVISION OF RECOMMENDATION ITU-R P.676-5 (29/11/2004)
%
% Created by Xiong Zubiao, 02/04/2005
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
deltao=(a5+a6.*theta).*(1e-4).*(p+e).*theta.^(0.8);
deltaw=0;
deltafo=a3.*1e-4.*(p.*theta.^(0.8-a4)+1.1.*e.*theta);
deltafw=b3.*1e-4.*(p.*theta.^b4+b5.*e.*theta.^b6);
%
%deltafo=(deltafo.^2+2.25e-6).^0.5;
%deltafw=0.535.*deltafw+(0.217.*deltafw.^2+2.1316e-12.*fw.^2./theta).^0.5;
%
part1Fio=(deltafo-deltao.*(fo-f))./((fo-f).^2+deltafo.^2);
part2Fio=(deltafo-deltao.*(fo+f))./((fo+f).^2+deltafo.^2);
%
Fio=(f./fo).*(part1Fio+part2Fio);
%
part1Fiw=(deltafw-deltaw.*(fw-f))./((fw-f).^2+deltafw.^2);
part2Fiw=(deltafw-deltaw.*(fw+f))./((fw+f).^2+deltafw.^2);
%
Fiw=(f./fw).*(part1Fiw+part2Fiw);
%
%
% Dry air continuum N2Df ( N''D(f) )
%
d=5.6e-4.*p.*theta.^0.8;
N2Df=f.*p.*theta.^2.*(6.14e-5./(d.*(1+(f./d).^2))+...
		(1.4e-12.*p.*theta.^1.5)./(1+1.9e-5.*f.^1.5));
%
% N2Wf=f.*(3.57.*theta.^7.5.*e+0.113.*p).*1.0e-7.*e.*theta.^3;
%
N2f=sum(Fio.*Sio)+sum(Fiw.*Siw)+N2Df;
%
y=0.1820.*f.*N2f;




