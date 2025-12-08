%function [e,TB,th]=RTModle_6(f,incident_angle,epsilon,mu,d,Tg,polarstyle)
function [eh,ev] = RefFourConceal(incidence,Param)
% 利用等效多层传输线模型，计算四层介质隐身目标的发射率
% 输入：入射角incidence
% 输出：输出发射率的水平极化eh和垂直极化分量ev
%%
Frequency = Param.Observe.center_freq;                 % 辐射计中心频率
Wavelength = 0.3./Frequency;                           % 波长，单位m
mu = zeros(1,5);         %各层%磁导率，不包含空气层
epsilon = zeros(1,5);     %各层介电常数，不包含空气层
thickness = zeros(1,4);   %各层厚度
%磁导率
mu(1)=Param.Target.layer1_target_permeability; 
mu(2) = Param.Target.layer2_target_permeability;          % 涂层1介电常数，红外2.3-0.02i
mu(3) = Param.Target.layer3_target_permeability;          % 涂层2介电常数，雷达20.6-3.07i
mu(4) = Param.Target.layer4_target_permeability;          % 涂层3介电常数，雷达18.4-2.57i
mu(5) = 1;
%介电常数
epsilon(1) = Param.Target.layer1_target_permittivity;          % 涂层1磁导率，红外0.67-0.75i
epsilon(2) = Param.Target.layer2_target_permittivity;          % 涂层2磁导率，雷达1.24-0.93i
epsilon(3)  = Param.Target.layer3_target_permittivity;          % 涂层3磁导率，雷达1.03-0.72i
epsilon(4)  =Param.Target.layer4_target_permittivity; 
epsilon(5) = -inf*1i;
%厚度
thickness(1) = Param.Target.layer1_coat_thickness;             % 涂层1厚度（单位m），红外0.000033
thickness(2)= Param.Target.layer2_coat_thickness;             % 涂层2厚度（单位m），雷达0.0005
thickness(3)= Param.Target.layer3_coat_thickness;             % 涂层3厚度（单位m），雷达0.0004
thickness(4)= Param.Target.layer4_coat_thickness;             % 涂层3厚度（单位m），雷达0.0004
%%
k0=2*pi/Wavelength;
k=zeros(1,5);  %存放各层 等效传播常数；
x=zeros(1,5);  %存放各层等效传播常数的实部
y=zeros(1,5);  %存放各层等效传播常数的虚部取反的值
p=zeros(1,5);  %存放各层等效传播常数的实部和虚部的积
kx=zeros(1,5);
q=zeros(1,5);
th=zeros(1,5);    %存储2-6的入射或反射角度
%% 
%计算各层的 等效传播常数 相关等参数
for i=1:5
 k(i)=k0*sqrt(mu(i)*epsilon(i));    %等效传播常数
 x(i)=real(k(i));                   %实部
 y(i)=-imag(k(i));                  %虚部
 p(i)=2*x(i)*y(i);                  %实部和虚部之积
end
%%
kx(1)=k0*sind(incidence);
q(1)=x(1)^2-y(1)^2-kx(1)^2;
th(1)=atand(kx(1)/(2^(-0.5)*((p(1)^2+q(1)^2)^(0.5)+q(1))^(0.5)));
for i=2:5
    kx(i)=k(i-1)*sind(th(i-1));
    q(i)=x(i)^2-y(i)^2-kx(i)^2;
    th(i)=real(atand(kx(i)/(2^(-0.5)*((p(i)^2+q(i)^2)^(0.5)+q(i))^(0.5))));
end
%%
eta=zeros(1,5);       %各层波阻抗
z=zeros(1,5);         %各层媒质的特征阻抗
gamma=zeros(1,5);
r=zeros(1,5);         %媒质i入射到i层边界处的反射系数
Zin=zeros(1,4);       %各层等效输入阻抗
%% 水平极化
%计算水平极化下各层波阻抗
for i=1:5
    eta(i)=sqrt(mu(i)./epsilon(i));%各层波阻抗
    z(i)=eta(i).*secd(th(i));%各层媒质的特征阻抗
    z(5)=-inf;
     gamma(i)=2*pi/Wavelength*sqrt(epsilon(i)); 
end
%计算水平极化媒质i入射到i层边界处的反射系数
for i=2:5
    r(i-1)=-(z(i)-z(i-1))./(z(i)+z(i-1));
end
r(4)=-1;
%由里层向外层逐步计算等效输入阻抗和反射系数
for i=4:-1:2
    Zin(i)=z(i)*((1+r(i)*exp(-2i.*gamma(i).*secd(th(i)).*thickness(i)))/(1-r(i)*exp(-2i.*gamma(i).*secd(th(i)).*thickness(i))));
    r(i-1)=(Zin(i)-z(i-1))./(Zin(i)+z(i-1));% 水平极化    
end
Zin(1)=z(1)*((1+r(1)*exp(-2i.*gamma(1).*secd(th(1)).*thickness(1)))./(1-r(1)*exp(-2i.*gamma(1).*secd(th(1)).*thickness(1))));
Rh=(Zin(1)-secd(incidence))./(Zin(1)+secd(incidence));%
eh=1-abs(Rh).^2;
%% 垂直极化
%计算垂直极化下各层波阻抗
for i=1:5
    eta(i)=sqrt(mu(i)./epsilon(i));%各层波阻
    z(i)=eta(i).*cosd(th(i));
    z(5)=-inf;
     gamma(i)=2*pi/Wavelength*sqrt(epsilon(i)); 
end
for i=2:5
    r(i-1)=-(z(i)-z(i-1))./(z(i)+z(i-1));
end
r(4)=-1;
%垂直极化时，由里层向外层逐步计算等效输入阻抗和反射系数
for i=4:-1:2
    Zin(i)=z(i)*((1+r(i)*exp(-2i.*gamma(i).*secd(th(i)).*thickness(i)))/(1-r(i)*exp(-2i.*gamma(i).*secd(th(i)).*thickness(i))));
    r(i-1)=-(Zin(i)-z(i-1))./(Zin(i)+z(i-1));
end
Zin(1)=z(1)*((1+r(1)*exp(-2i.*gamma(1).*secd(th(1)).*thickness(1)))./(1-r(1)*exp(-2i.*gamma(1).*secd(th(1)).*thickness(1))));
Rv=-(Zin(1)-cosd(incidence))./(Zin(1)+cosd(incidence));%
ev=1-abs(Rv).^2;
end
