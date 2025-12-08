function Kl = GetCoefficient (T,f)
%计算云、雾条件下的比衰减系数
%参考ITU-R P.840-5

theta = 300/T;
permit0 = 77.6 + 103.3*(theta-1);
permit1 = 5.48;
permit2 = 3.51;

fp = 20.09 - 142*(theta-1) + 294*(theta-1)^2;
fs = 590 - 1500*(theta -1);

permit_2 = f*(permit0-permit1)/fp/(1+(f/fp)^2) + f*(permit1-permit2)/fs/(1+(f/fs)^2);
permit_1 = (permit0-permit1)/(1+(f/fp)^2) + (permit1-permit2)/(1+(f/fs)^2) + permit2;
eta = (permit_1+2)/permit_2;
Kl = 0.819*f/permit_2/(1+eta^2);
end