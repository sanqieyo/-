function y = GetWeatherAtt(weather,T,f,height)
%计算特定云或雾条件下的具体衰减量
%参考ITU-R P.840-5与ITU-R P.838-2
%ro为云或雾中的液体水密度,g/m3
%rainrate降雨率
y=0;

%云
if strcmp(weather,'CLOUDY')
    if (height>=1)&&(height<=4)
        ro=0.1;
        y=GetCoefficient(f,T)*ro;
    end
end

%雾
if strcmp(weather,'MIST')
    if  height<=0.1
        ro=0.156;
        y=GetCoefficient(f,T)*ro;
    end
end

%雨
if strcmp(weather,'RAINY')
    if  height<=5
        rainrate=1;
        tao=45;
        y=getrainattunatione(f,rainrate,tao);
    end
end

end