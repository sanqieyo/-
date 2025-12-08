function sea_point = oceanprodct(Param,tt)
%%%%%%%二维海面产生%%%%%%%%%%%

% 本函数用来根据海浪谱生成海面及海面的斜率

% 输入变量：
% x_length     海面x轴方向尺寸，单位为米
% y_length     海面y轴方向尺寸，单位为米
% dx           海面面元尺寸，单位为米
% dy           海面面元尺寸，单位为米
% t 时间
% ocean_type   海浪谱类型：
%              ocean_type = 'PM', 表示PM谱
%              ocean_type = 'JONS', 表示JONS谱
%              ocean_type = 'DV', 表示DV谱
%              ocean_type = 'EL', 表示EL谱
%fine_D        全局观测角

% 输出变量：
% height         海面位移
% slope_x        海面横坐标方向斜率
% slope_y        海面纵坐标方向斜率

%--------------------初始化------------------------------------%
x_length = Param.x_length;
y_length = Param.y_length;
dx = Param.delta_x; 
dy = Param.delta_y;
t = tt;
ocean_type = 'EL';
fine_D = 0;
fine_wind = Param.wind_angle;
wind_speed = Param.wind_speed;


%-----------------------------------------Start 设置变量-----------------------------------------%
%-------------------剖分网格格子横纵坐标数目-------------------%
M=ceil(x_length/dx);
N=ceil(y_length/dy);

%-------------------波数空间的分辨率-------------------%
dkx=2*pi/x_length;
dky=2*pi/y_length;

%-------------------离散网格下的离散坐标-------------------%
m=0:M;%波数横坐标轴序列
n=0:N;%波数纵坐标轴序列
p=0:M;%空间横坐标轴序列
q=0:N;%空间纵坐标轴序列

mid=ones(1,N+1);%中间变量
m_=m'*mid;%波数空间中的波数横坐标轴序列值
n_=mid'*n;%波数空间中的波数纵坐标轴序列值
p_=m_;%空间中的空间横坐标轴序列值
q_=n_;%空间中的空间纵坐标轴序列值

%-------------------离散坐标对应的离散波数值-------------------%
kx=-pi/dx+m*2*pi/x_length;%离散波数空间横坐标
ky=-pi/dy+n*2*pi/y_length;%离散波数空间纵坐标
    kx_odds=kx'*mid;%组合排列矩阵
    ky_odds=mid'*ky;%组合排列矩阵
    kx_odds1=kx_odds;
    ky_odds1=ky_odds;
    k_combined=(kx_odds.^2+ky_odds.^2).^0.5;%合成波数矢量的大小
x=p*x_length/M-x_length/2;%离散空间横坐标
y=q*y_length/N-y_length/2;%离散空间纵坐标
    x_odds=x'*mid;%组合排列矩阵
    y_odds=mid'*y;%组合排列矩阵

%-------------------随机初相位以及传播方位生成-------------------%
% 海况条件设置(请根据给定的海况类型ocean_type设置相应)
gravity=9.8;%海面重力加速度
delta_phi=zeros(M+1,N+1);
delta_phi(1:M/2,N/2+1:N+1)=pi;
delta_phi(1:M/2,1:N/2)=pi;
delta_phi(M/2+1:M+1,1:N/2)=2*pi;
fine_transmit=delta_phi+atan((ky_odds./kx_odds));%传播方位变量

random_sigma = randomMN(M,N);%随机相位生成
wn=(k_combined*gravity).^0.5;%角频率
%-----------------------------------------END 设置变量-----------------------------------------%

%--------------------------------Start 生成波高谱与归一化振幅谱-----------------------------------%
%方法与公式参考专利“基于蒙特卡罗法的单向传播的海浪海面生成方法”
L = length(ocean_type);
F = 0;
amn = 0;
for H=1:L
    [F_,fine_thistype] = hei_spec(kx_odds,ky_odds,k_combined, fine_D, ocean_type, 'F2D',wind_speed,fine_wind);
    amn_=2*sqrt(F_.*dkx.*dky);
    amn_(isnan(amn_))=0;
    [ amn_,F_ ] = zerolize( amn_,F_,fine_D,fine_thistype,kx_odds,ky_odds );
    amn = amn+amn_;
    F = F+F_;
end
%--------------------------------End 生成波高谱与归一化振幅谱-----------------------------------%

%---------------------------Start 生成海面与海面斜率---------------------------%
posneg=exp(1i*pi*(m_+n_));
middle=random_sigma.*amn.*posneg.*exp(1i.*(-wn*t));
posneg1=posneg(2:end,2:end);
middle1=middle(2:end,2:end);
middle_sx=middle1.*kx_odds(2:end,2:end)*1i;
middle_sy=middle1.*ky_odds(2:end,2:end)*1i;
height=M*N*ifft2(middle1).*posneg1;
slope_x=M*N*ifft2(middle_sx).*posneg1;
slope_y=M*N*ifft2(middle_sy).*posneg1;
sea_point=real(height);
slope_x=real(slope_x);
slope_y=real(slope_y);
%-----------------------------------------END 生成海面与海面斜率-----------------------------------------%
end




