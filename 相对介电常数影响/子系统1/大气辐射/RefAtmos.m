function y=RefAtmos(Ph,T,ro)
%
% FUNCTION Y=REF_IND(Ph,T,ro)
%
% Non-dispersive Refractivy (frequency independemt component of the
% refractivity in ppm (parts per million)
%
% Ph = Pressure in hPa
% p= dry air pressure in KPa
% e= water-vapour partial pressure in Kpa
% P=total pressure = p+e in Kpa
% T= temperature in Kelvin
% ro= water vapour concentration (g/m^3)
%
% 24/6/95 OK
%
% J.P.V. Poiares Baptista
% ESA/ESTEC, Keplerlaan 1
% NL-2200 AG Noordwijk
%
P=Ph./10.0;
theta=300.0./T;
%e=ro./(7.223.*theta);
e=ro.*T./216.7;
p=P-e;
%
% Dry air component
%
Nd=2.588.*p.*theta;
%
% Water vapour component
%
Nw=(41.63.*theta+2.39).*e.*theta;
%
y=Nd+Nw;


