function [eh,ev] = RefTwoConceal(incidence,Param)
% 利用等效多层传输线模型，计算双层介质隐身目标的发射率
% 输入：入射角incidence
% 输出：输出发射率的水平极化eh和垂直极化分量ev

%% 涂层参数

Frequency = Param.Observe.center_freq;                 % 辐射计中心频率
Wavelength = 0.3./Frequency;                           % 波长，单位m
jd1=1;
jd2 = Param.Target.layer1_target_permittivity;          % 涂层1介电常数，红外2.3-0.02i
jd3 = Param.Target.layer2_target_permittivity;          % 涂层2介电常数，雷达20.6-3.07i
jd4=-inf.*1i;

cd1=1;
cd2 = Param.Target.layer1_target_permeability;          % 涂层1磁导率，红外0.67-0.75i
cd3 = Param.Target.layer2_target_permeability;          % 涂层2磁导率，雷达1.24-0.93i
cd4=1;

thickness1 = Param.Target.layer1_coat_thickness;             % 涂层1厚度（单位m），红外0.000033
thickness2 = Param.Target.layer2_coat_thickness;             % 涂层2厚度（单位m），雷达0.0005

incidence_2=acosd(sqrt(1-((sind(incidence)).^2)./(cd2.*jd2)));
incidence_3=acosd(sqrt(1-((sind(incidence)).^2)./(cd3.*jd3)));
incidence_4=acosd(sqrt(1-((sind(incidence)).^2)./(cd4.*jd4)));

%% 水平极化p=h,n=nh=0

nh=0;
Z1=sqrt(cd1./jd1).*secd(incidence);
Z2=sqrt(cd2./jd2).*secd(incidence_2);
Z3=sqrt(cd3./jd3).*secd(incidence_3);
Z4=sqrt(cd4./jd4).*secd(incidence_4);
R1=((-1).^nh).*((Z2-Z1)./(Z2+Z1));
R2=((-1).^nh).*((Z3-Z2)./(Z3+Z2));
R3=((-1).^nh).*((Z4-Z3)./(Z4+Z3));
r2=(2.*pi./Wavelength).*(sqrt(cd2.*jd2));
r3=(2.*pi./Wavelength).*(sqrt(cd3.*jd3));

Re2=(R2+R3.*exp(-2i.*thickness2.*r3.*secd(incidence_3)))./(1+R2.*R3.*exp(-2i.*thickness2.*r3.*secd(incidence_3)));   % 化简式
Reh=(R1+Re2.*exp(-2i.*thickness1.*r2.*secd(incidence_2)))./(1+R1.*Re2.*exp(-2i.*thickness1.*r2.*secd(incidence_2))); % 化简式
eh=1-abs(Reh).^2;               % 发射率水平极化分量

%% 垂直极化p=v,n=nv=1

nv=1;
Z1=sqrt(cd1./jd1).*cosd(incidence);
Z2=sqrt(cd2./jd2).*cosd(incidence_2);
Z3=sqrt(cd3./jd3).*cosd(incidence_3);
Z4=sqrt(cd4./jd4).*cosd(incidence_4);
R1=((-1).^nv).*((Z2-Z1)./(Z2+Z1));
R2=((-1).^nv).*((Z3-Z2)./(Z3+Z2));
R3=((-1).^nv).*((Z4-Z3)./(Z4+Z3));
r2=(2.*pi./Wavelength).*(sqrt(cd2.*jd2));
r3=(2.*pi./Wavelength).*(sqrt(cd3.*jd3));

Re2=(R2+R3.*exp(-2i.*thickness2.*r3.*secd(incidence_3)))./(1+R2.*R3.*exp(-2i.*thickness2.*r3.*secd(incidence_3)));   % 化简式
Rev=(R1+Re2.*exp(-2i.*thickness1.*r2.*secd(incidence_2)))./(1+R1.*Re2.*exp(-2i.*thickness1.*r2.*secd(incidence_2))); % 化简式
ev = 1-abs(Rev).^2;           % 发射率垂直极化分量

end