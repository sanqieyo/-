function [p,T,ro] = StdAtmos(atmos,hei,ground_temp)
%
%function y=[p, T ro]=st_atmos(atmos,hei)
%
% FUNCTION TO DERIVE THE PROFILES OF PRESSURE (p in hPa), TEMPERATURE (T in K)
% AND WATER VAPOUR CONCENTRATION (ro in g/m^3)
%
% y=[p,T,ro]
% atmos='SUMMER MIDLAT','WINTER MIDLAT','SUMMER SUBART',  
%			'WINTER SUBART' OR 'ANNUAL TROPIC',
%			'ITU-R REF STD'
% hei= height in kilometers   ¸ß¶È£¨km£©
%
% CAREFUL: use only with one value for the height at a time!!!
%
% bug detected in the ITU standard atmosphere
% corrected 18/9/96 - Pedro
%
% corrected wrong sign detected in SUMMER SUBART atmosphere
% 12/1/97 - Pedro
%
%  5/3/95 - now includes also ITU-R PN Rec. 835
%
%  28/9/96 - found another bug in ITU atmosphere - pressure corrected
%
% J.P.V. Poiares Baptista
% ESA/ESTEC, Keplerlaan 1
% NL-2200 AG Noordwijk
%
if strcmpi(atmos,'SUMMER MIDLAT')
	% Pressure in hPa
	if hei>=0 && hei<=10,
		p=1012.8186-111.5569.*hei+3.8646.*hei.^2;
	elseif hei>10 && hei<=72,
		p=283.7096.*exp(-0.147.*(hei-10));
	elseif hei>72 && hei<=100,
		p=0.0312402229.*exp(-0.165.*(hei-72));
	else
		p=0;
	end
	% Temperature in Kelvin
	if hei>=0 && hei<=13,
		T=ground_temp-5.2159.*hei-0.07109.*hei.^2;          % 294.9838
	elseif hei>13 && hei<=17,
		T=215.15;
	elseif hei>17 && hei<=47,
		T=215.15.*exp(0.008128.*(hei-17));
	elseif hei>47 && hei<=53,
		T=275;
	elseif hei>53 && hei<=80,
		T=275+20.*(1-exp(0.06.*(hei-53)));
	elseif hei>80 && hei<=100,
		T=175;
	else
		T=0;
	end
	% absolute humidity in g/m^3
	if hei>=0 && hei<=15,
		ro=14.2542.*exp(-0.4174.*hei-0.02290.*hei.^2+0.001007.*hei.^3);
	else
		ro=0;
	end
	%
elseif strcmpi(atmos,'WINTER MIDLAT')
	% Pressure in hPa
	if hei>=0 && hei<=10,
		p=1018.8627-124.2954.*hei+4.8307.*hei.^2;
	elseif hei>10 && hei<=72,
		p=258.9787.*exp(-0.147.*(hei-10));
	elseif hei>72 && hei<=100,
		p=0.0285170199.*exp(-0.155.*(hei-72));
	else
		p=0;
	end
	% Temperature in Kelvin
	if hei>=0 && hei<=10,
		T=272.7241-3.6217.*hei-0.1759.*hei.^2;
	elseif hei>10 && hei<=33,
		T=218;
	elseif hei>33 && hei<=47,
		T=218+3.3571.*(hei-33);
	elseif hei>47 && hei<=53,
		T=265;
	elseif hei>53 && hei<=80,
		T=265-2.0370.*(hei-53);
	elseif hei>80 && hei<=100,
		T=210;
	else
		T=0;
	end
	% absolute humidity in g/m^3
	if hei>=0 && hei<=10,
		ro=3.4742.*exp(-0.2697.*hei-0.03604.*hei.^2+0.0004489.*hei.^3);
	else
		ro=0;
	end
	%
elseif strcmpi(atmos,'SUMMER SUBART')
	% Pressure in hPa
	if hei>=0 && hei<=10,
		p=1008.0278-113.2494.*hei+3.9408.*hei.^2;
	elseif hei>10 && hei<=72,
		p=269.6138.*exp(-0.140.*(hei-10));
	elseif hei>72 && hei<=100,
		p=0.0458211532.*exp(-0.165.*(hei-72));
	else
		p=0;
	end
	% Temperature in Kelvin
	if hei>=0 && hei<=10,
		T=286.8374-4.7805.*hei-0.1402.*hei.^2;
	elseif hei>10 && hei<=23,
		T=225;
	elseif hei>23 && hei<=48,
		T=225.*exp(0.008317.*(hei-23));
	elseif hei>48 && hei<=53,
		T=277;
	elseif hei>53 && hei<=79,
		T=277-4.0769.*(hei-53);
	elseif hei>79 && hei<=100,
		T=171;
	else
		T=0;
	end
	% absolute humidity in g/m^3
	if hei>=0 && hei<=15,
		ro=8.988.*exp(-0.3614.*hei-0.005402.*hei.^2-0.001955.*hei.^3);
	else
		ro=0;
	end
	%
elseif strcmpi(atmos,'WINTER SUBART')
	% Pressure in hPa
	if hei>=0 &&	hei<=10,
		p=1010.8828-122.2411.*hei+4.554.*hei.^2;
	elseif hei>10 && hei<=72,
		p=243.8718.*exp(-0.147.*(hei-10));
	elseif hei>72 && hei<=100,
		p=0.0268535481.*exp(-0.150.*(hei-72));
	else
		p=0;
	end
	% Temperature in Kelvin
	if hei>=0 && hei<=8.5,
		T=257.4345+2.3474.*hei-1.5479.*hei.^2+0.08473.*hei.^3;
	elseif hei>8.5 && hei<=30,
		T=217.5;
	elseif hei>30 && hei<=50,
		T=217.5+2.125.*(hei-30);
	elseif hei>50 && hei<=54,
		T=260;
	elseif hei>54 && hei<=100,
		T=260-1.667.*(hei-54);
	else
		T=0;
	end
	% absolute humidity in g/m^3
	if hei>=0 && hei<=10,
		ro=1.2319.*exp(0.07481.*hei-0.0981.*hei.^2+0.00281.*hei.^3);
	else
		ro=0;
	end
	%
elseif strcmpi(atmos,'ANNUAL TROPIC')
	% Pressure in hPa
	if hei>=0 &&hei<=10,
		p=1012.0306-109.0338.*hei+3.6316.*hei.^2;
	elseif hei>10 && hei<=72,
		p=284.8526.*exp(-0.147.*(hei-10));
	elseif hei>72 && hei<=100,
		p=0.0313660825.*exp(-0.165.*(hei-72));
	else
		p=0;
	end
	% Temperature in Kelvin
	if hei>=0 && hei<=17,
		T=300.4222-6.3533.*hei+0.005886.*hei.^2;
	elseif hei>17 && hei<=47,
		T=194+2.533.*(hei-17);
	elseif hei>47 && hei<=52,
		T=270;
	elseif hei>52 && hei<=80,
		T=270-3.0714.*(hei-52);
	elseif hei>80 && hei<=100,
		T=184;
	else
		T=0;
	end
	% absolute humidity in g/m^3
	if hei>=0 && hei<=15,
		ro=19.6542.*exp(-0.2313.*hei-0.1122.*hei.^2+0.01351.*hei.^3.0...
			-0.0005923.*hei.^4);
	else
		ro=0;
	end
	%
elseif strcmpi(atmos,'ITU-R REF STD')
	% Pressure in hPa
	if hei>=0 && hei<=11,
		p=1013.25.*(288.15./(288.15+(-6.5).*(hei))).^(34.163./(-6.5));
	elseif hei>11 && hei<=20,
		p=226.3226.*exp((-34.163.*(hei-11))./216.65);
	elseif hei>20 && hei<=32,
		p=54.7498.*(216.65./(216.65+(1.0).*(hei-20))).^(34.163);
	elseif hei>32 && hei <=47,
		p=8.6804.*(228.65./(228.65+(2.8).*(hei-32))).^(34.163./(2.8));
	elseif hei>47 && hei<=51,
		p=1.1091.*exp((-34.163.*(hei-47))./270.65);
	elseif hei>51 && hei<=71,
		p=0.6694.*(270.65./(270.65+(-2.8).*(hei-51))).^(34.163./(-2.8));
	elseif hei>71 && hei<=100,
		p=0.0396.*(214.65./(214.65+(-2.0).*(hei-71))).^(34.163./(-2.0));
	else
		p=0;
	end
	% Temperature in Kelvin
	if hei>=0 && hei<=11,
		T=288.15+(-6.5).*(hei);
	elseif hei>11 && hei<=20,
		T=216.65;
	elseif hei>20 && hei<=32,
		T=216.65+(hei-20);
	elseif hei>32 && hei <=47,
		T=228.65+2.8.*(hei-32);
	elseif hei>47 && hei<=51,
		T=270.65;
	elseif hei>51 && hei<=71,
		T=270.65-2.8.*(hei-51);
	elseif hei>71 && hei<=100,
		T=214.65-2.0.*(hei-71);
	else
		T=0;
	end
	% absolute humidity in g/m^3
	if hei>100,
		ro=0;
	else
		ro=7.5.*exp(-0.5.*hei);
		mix=(ro.*T./216.7)./p;
		if mix<2.0e-6,
			e=2.0e-6.*p;
			ro=e.*216.7./T;
		end
		
	end
	%
else
end



