function [yis_ZV, yis_ZVD, yis_EI] = Input_shape(ts,y1,Fs,b)
Ts = ts(2)-ts(1);

wn = 2*pi*Fs; % w=2*pi*f, 系统频率
wd = wn*sqrt(1-b*b); %阻尼振荡频率
T_Z = pi/wd; % 峰值时间
T_EI = pi/wn;
% 影响模型精度的参数：
% 输入分离位置tN,和阻尼系数b
f = 1/Ts; 
tN_Z = fix(f*T_Z);% 取整函数
tN_EI = fix(f*T_EI);

y2_Z = [zeros(tN_Z,1);y1(1:end-tN_Z,1)];
y3_Z = [zeros(2*tN_Z,1);y1(1:end-2*tN_Z,1)];
y2_EI = [zeros(tN_EI,1);y1(1:end-tN_EI,1)];
y3_EI = [zeros(2*tN_EI,1);y1(1:end-2*tN_EI,1)];
yis_ZV = [];
yis_ZVD = [];
yis_EI = [];
for n = 1:length(y1)
    yis_ZV(n,1) = fcnZV(y1(n),y2_Z(n),b);
    yis_ZVD(n,1) = fcnZVD(y1(n),y2_Z(n),y3_Z(n),b);
    yis_EI(n,1) = fcnEI(y1(n),y2_EI(n),y3_EI(n));
end
figure;
stairs(ts,y1,'color','m','LineWidth',1.2);hold on;
stairs(ts,yis_ZV,'-.','color','r','LineWidth',1.2);hold on;
stairs(ts,yis_ZVD,'-.','color','g','LineWidth',1.2);hold on;
stairs(ts,yis_EI,'-.','color','b','LineWidth',1.2);hold on;
legend('原输入','ZV整形曲线','ZVD整形曲线','EI整形曲线')
xlabel('Time(s)');ylabel('Amp')
grid on;
end

function y = fcnEI(u1,u2,u3)
    Vexp = 0.01; %百分比
    A1 = (1+Vexp)/4;
    A2 = (1-Vexp)/2;
    A3 = (1+Vexp)/4;
    y = A1*u1 + A2*u2 + A3*u3;
end

function y = fcnZV(u1,u2,bs)
%     bs = 0.1; 
    K = exp(-bs*pi/(sqrt(1-bs*bs))); 
    
    A1 = 1/(1+K);
    A2 = K/(1+K);
    y = A1*u1 + A2*u2;
end

function y = fcnZVD(u1,u2,u3,bs)
%     bs = 0.1; 
    K = exp(-bs*pi/(sqrt(1-bs*bs))); 
    
    A1 = 1/(1+K)^2;
    A2 = 2*K/(1+K)^2;
    A3 = K^2/(1+K)^2;
    y = A1*u1 + A2*u2 + A3*u3;
end