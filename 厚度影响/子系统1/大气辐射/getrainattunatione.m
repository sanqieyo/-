function y=getrainattunatione(f,rainrate,tao)
% getrainattunatione 求雨的衰减系数，单位dB/km。
% f为频率，单位GHz；rainrate为降雨率，单位mm/h；tao为极化率，如45度为圆极化；sita为天线仰角，单位度。
% clear;
% rainrate=0.25;
% tao=45;
% sita=45;%测试
global SimuParam
sita = SimuParam.elevation;
sita=sita./180.*pi;
tao=tao./180.*pi;
kha=[-5.33980 -0.35351 -0.23789 -0.94158];
khb=[-0.10008 1.26970 0.86036 0.64552];
khc=[1.13098 0.45400 0.15354 0.16817];
khmk=-0.18961;
khck=0.71147;

kva=[-3.80595 -3.44965 -0.39902 0.50167];
kvb=[0.56934 -0.22911 0.73042 1.07319];
kvc=[0.81061 0.51059 0.11899 0.27195];
kvmk=-0.16398;
kvck=0.63297;

alphaha=[-0.14318 0.29591 0.32177 -5.37610 16.1721];
alphahb=[1.82442 0.77564 0.63773 -0.96230 -3.29980];
alphahc=[-0.55187 0.19822 0.13164 1.47828 3.43990];
alphahma=0.67849;
alphahca=-1.95537;

alphava=[-0.07771 0.56727 -0.20238 -48.2991 48.5833];
alphavb=[2.33840 0.95545 1.14520 0.791669 0.791459];
alphavc=[-0.76284 0.54039 0.26809 0.116226 0.116479];
alphavma=-0.053739;
alphavca=0.83433;
% f=1:1000;
n=length(f);
for i=1:n
    kh(i)=10.^(sum(kha.*exp(-((log10(f(i))-khb)./khc).^2))+khmk.*log10(f(i))+khck);
    alphah(i)=sum(alphaha.*exp(-((log10(f(i))-alphahb)./alphahc).^2))+alphahma.*log10(f(i))+alphahca;

    kv(i)=10.^(sum(kva.*exp(-((log10(f(i))-kvb)./kvc).^2))+kvmk.*log10(f(i))+kvck);
    alphav(i)=sum(alphava.*exp(-((log10(f(i))-alphavb)./alphavc).^2))+alphavma.*log10(f(i))+alphavca;
end

k=(kh+kv+(kh-kv).*(cos(sita)).^2.*cos(2.*tao))./2;
alpha=(kh.*alphah+kv.*alphav+(kh.*alphah-kv.*alphav).*(cos(sita)).^2.*cos(2.*tao))./2./k;

rainattunatione=k.*rainrate.^alpha;
y=rainattunatione;
end