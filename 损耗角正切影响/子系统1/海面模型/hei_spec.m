function [ FD,fine_thistype] = hei_spec( kx_odds,ky_odds,K, fine_D, ocean_type, spec_type, wind_speed,fine_wind )
% spec_type的可取范围为[F1D,S1D,F2D]，分别表示某个方位向的波高谱、全向谱、二维波高谱
% K的维度应当与spec_type指定的维度一样。例如spec_type=F1D,则K为一维。
% F1D是一维波高谱，S1D是一维全向谱，F2D是二维波高谱
[M_,N_]=size(K);
M=M_-1;
N=N_-1;
% 请在这里更改海况参数设置
% fine_wind = 0*pi/180;%全局风向角
fine_wind = fine_wind*pi/180;
% wind_speed = 12;%风速
omega_ = 0.84;%风浪逆波龄

Hs_ = 4;%涌浪有效波高
kpdv = 2*pi/200;%涌浪峰值波数
fine_DV = 5*pi/180;%涌浪峰值波向


gravity=9.81;
if strcmp(ocean_type,'PM')
%---------------------------Start PM海浪模型波高谱---------------------------%
    k_peak=0.7694*gravity/(wind_speed^2);
    fine_thistype = fine_wind;
    if strcmp(spec_type(end-1:end),'1D')
        % 定义波高谱变量为：
        FD=exp(-5/4*(k_peak./K).^2)*(0.008/2).*K.^(-3);%PM谱波高谱变量
        if strcmp(spec_type(1),'F')
            FD = FD ./ K .* 4.*cos(fine_thistype-fine_D).^4./(3*pi);
        end
    elseif strcmp(spec_type(end-1:end),'2D')
        delta_phi=zeros(M+1,N+1);
        delta_phi(1:M/2,N/2+1:N+1)=pi;
        delta_phi(1:M/2,1:N/2)=pi;
        delta_phi(M/2+1:M+1,1:N/2)=2*pi;

        fine_transmit=delta_phi+atan((ky_odds./kx_odds));
        FD=exp(-5/4*(k_peak./K).^2)*(0.0081/2).*K.^(-4).*4.*cos((fine_thistype-fine_D)-fine_transmit).^4./(3*pi);%PM谱波高谱变量
    end
%---------------------------END PM海浪模型波高谱---------------------------%

elseif strcmp(ocean_type,'JONS')
%---------------------------Start JONS海浪模型波高谱---------------------------%
    kafang=90000;%60km风区
    kafang_w=gravity*kafang/wind_speed^2;
    alfa=0.076*(kafang_w)^(-0.22);
    gama=3.3;
    w0_w=18.3*kafang_w^(-1/3);
    w0=w0_w*gravity/wind_speed;
    k_peak=w0^2/gravity;
    fine_thistype = fine_wind;
    if strcmp(spec_type(end-1:end),'1D')
        deltas_l=zeros(1,N+1);
        deltas_s=zeros(1,N+1);
        [row1,col1]=find(K>k_peak);
        for hh1=1:length(row1)
            deltas_l(row1(hh1),col1(hh1))=0.09;
        end
        
        [row2,col2]=find(K<=k_peak);
        for hh1=1:length(row2)
            deltas_s(row2(hh1),col2(hh1))=0.07;
        end
        
        deltas=deltas_l+deltas_s;
        FD=alfa/2.*K.^(-3).*exp(-5/4.*(k_peak./K).^2).*gama.^exp(-((K./k_peak).^0.5-1).^2./(2.*deltas.^2));
        if strcmp(spec_type(1),'F')
            FD = FD ./ K .* 4.*cos(fine_thistype-fine_D).^4./(3*pi);
        end
    elseif strcmp(spec_type(end-1:end),'2D')
        deltas_l1=zeros(M+1,N+1);
        deltas_s1=zeros(M+1,N+1);
        [row1,col1]=find(K>k_peak);
        for hh1=1:length(row1)
            deltas_l1(row1(hh1),col1(hh1))=0.09;
        end
        [row2,col2]=find(K<=k_peak);
        for hh1=1:length(row2)
            deltas_s1(row2(hh1),col2(hh1))=0.07;
        end
        deltas1=deltas_l1+deltas_s1;
        delta_phi=zeros(M+1,N+1);
        delta_phi(1:M/2,N/2+1:N+1)=pi;
        delta_phi(1:M/2,1:N/2)=pi;
        delta_phi(M/2+1:M+1,1:N/2)=2*pi;

        fine_transmit=delta_phi+atan((ky_odds./kx_odds));
        FD=alfa/2.*K.^(-4).*exp(-5/4.*(k_peak./K).^2).*gama.^exp(-((K./k_peak).^0.5-1).^2./(2.*deltas1.^2)).*4.*cos((fine_thistype-fine_D)-fine_transmit).^4./(3*pi);
    end
elseif strcmp(ocean_type,'DV')
%---------------------------START2D DV海浪模型波高谱---------------------------%
    fine_thistype = fine_DV;
    deltas=0.006;
    Hs=Hs_;%%%有效波高设置
    k_peak=kpdv;
    FD_=Hs.^2./(16*sqrt(2*pi)*deltas).*exp(-0.5.*((K-k_peak)/deltas).^2);
    
    
    syms xx yy;
    cp=cos(xx-yy).^14;
    cpp=int(cp,xx,0,2*pi);
    
    F322D1=cos((fine_thistype-fine_D)).^14./double(cpp);
    
    if strcmp(spec_type(end-1:end),'1D')
        FD = FD_;
        if strcmp(spec_type(1),'F')
            FD = FD ./ K .* F322D1;
        end
    elseif strcmp(spec_type(end-1:end),'2D')
        delta_phi=zeros(M+1,N+1);
        delta_phi(1:M/2,N/2+1:N+1)=pi;
        delta_phi(1:M/2,1:N/2)=pi;
        delta_phi(M/2+1:M+1,1:N/2)=2*pi;

        fine_transmit=delta_phi+atan((ky_odds./kx_odds));
        F322D=cos((fine_thistype-fine_D)-fine_transmit).^14./double(cpp);
        FD=FD_./K.*F322D;
    end
%---------------------------END2D DV海浪模型波高谱---------------------------%
elseif strcmp(ocean_type,'EL')
    fine_thistype = fine_wind;
    omega = omega_;
    km = 363;
    cm = sqrt(2*gravity/km); %HP师姐的，文献9给的是cm=0.23m/s
%---------------------------Start EL海浪模型波高谱---------------------------%
    %----------Start 计算高频曲率谱----------%
    ap = 0.006.*sqrt( omega );
    k0 = gravity ./ (wind_speed.^2);
    kp = k0 .* omega.^2;
    cp = wind_speed/omega;%来自文献9，公式32上一段文字
    sigma = 0.08*(1+4*omega.^(-3));
    [~,us,~,~,~]=blpara(wind_speed);
        
    if us<cm
        am = 0.01*(1+log(us/cm));
    elseif us>cm
        am = 0.01*(1+3*log(us/cm));
    end
    if omega>=0.84 && omega<1
        gama = 1.7;
    elseif omega>1 && omega<5
        gama = 1.7 + 6*log10(omega);
    end
%     cp = sqrt(gravity ./ kp);% 来自HP师姐的论文
    if strcmp(spec_type(end-1:end),'1D')
        cc = sqrt(gravity./K.*(1+(K/km).^2));
        Lpm = exp( -5/4 * (kp./K).^2 );
        
        
        great_gama = exp(-(sqrt(K./kp)-1).^2./(2*sigma.^2));
        
        Jp = gama.^great_gama;
        
        Bl = ap/2 * cp ./ cc .* Lpm .* Jp .* exp(-omega/sqrt(10)*(sqrt(K./kp)-1));
        %----------END 计算高频曲率谱----------%
        
        %----------Start 计算低频曲率谱----------%
        Bh = am./2.*cm./cc.*exp(-0.25*(K/km-1).^2);
        %----------END 计算低频曲率谱----------%
        
        Sk = K.^(-3).*(Bl+Lpm.*Bh);
        
        Fk = Sk./K;% 计算二维谱使用的Fk
        
        [deltak]=tony_spread(K,wind_speed,omega);
        %     F1D = Fk./2/pi.*(1+deltak.*cos(2.*(fine_wind-fine_watch)));
        FD = Fk.*K;
        if strcmp(spec_type(1),'F')
            FD = FD ./ K ./2/pi.*(1+deltak.*cos(2.*(fine_thistype-fine_D)));
        end
    elseif strcmp(spec_type(end-1:end),'2D')
        cc1 = sqrt(gravity./K.*(1+(K/km).^2));
  
        Lpm1 = exp( -5/4 * (kp./K).^2 );
        
        great_gama1 = exp(-(sqrt(K./kp)-1).^2./(2*sigma.^2));
        Jp1 = gama.^great_gama1;
        
        Bl1 = ap/2 * cp ./ cc1 .* Lpm1 .* Jp1 .* exp(-omega/sqrt(10)*(sqrt(K./kp)-1));
        %----------END 计算高频曲率谱----------%
        
        %----------Start 计算低频曲率谱----------%
        Bh1 = am/2*cm./cc1.*exp(-0.25*(K/km-1).^2);
        %----------END 计算低频曲率谱----------%
        
        %     Sk = K.^(-3).*(Bl + Bh);%全向谱Sk
        Sk1 = K.^(-3).*(Bl1+Lpm1.*Bh1);
        
        Fk1 = Sk1./K;% 计算二维谱使用的Fk
        
        [deltak1]=tony_spread(K,wind_speed,omega);
        delta_phi=zeros(M+1,N+1);
        delta_phi(1:M/2,N/2+1:N+1)=pi;
        delta_phi(1:M/2,1:N/2)=pi;
        delta_phi(M/2+1:M+1,1:N/2)=2*pi;

        fine_transmit=delta_phi+atan((ky_odds./kx_odds));
        
        FD = Fk1./2/pi.*(1+deltak1.*cos(2.*((fine_thistype-fine_D)-fine_transmit)));
    end
end
end