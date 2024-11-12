clc;clear;close all;

Ts = 1/4000;
[y1,ts] = gensig('sq',0.2,1,Ts);

Fs = 40; %Hz
wn = 2*pi*Fs; % w=2*pi*f, 系统频率
b = 0.1;
G = tf(wn^2,[1 2*b*wn wn^2]);
wd = wn*sqrt(1-b*b); %阻尼振荡频率
T_Z = pi/wd; % 峰值时间
T_EI = pi/wn;
% 影响模型精度的参数：
% 输入分离位置tN,和阻尼系数b
f = 4000; 
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
    yis_ZV(n,1) = fcnZV(y1(n),y2_Z(n));
    yis_ZVD(n,1) = fcnZVD(y1(n),y2_Z(n),y3_Z(n));
    yis_EI(n,1) = fcnEI(y1(n),y2_EI(n),y3_EI(n));
end
[input_ZV, input_ZVD, input_EI] = Input_shape(ts,y1,Fs,b);
figure;
stairs(ts,y1,'color','m','LineWidth',1.2);hold on;
stairs(ts,yis_ZV,'-.','color','r','LineWidth',1.2);hold on;
stairs(ts,yis_ZVD,'-.','color','g','LineWidth',1.2);hold on;
stairs(ts,yis_EI,'-.','color','b','LineWidth',1.2);hold on;
legend('方波','ZV整形曲线','ZVD整形曲线','EI整形曲线')
xlabel('Time(s)');ylabel('Amp')
ylim([0 1.2])
grid on;

[out, ts] = lsim(G,y1,ts);
[outZV, ts] = lsim(G,yis_ZV,ts);
[outZVD, ts] = lsim(G,yis_ZVD,ts);
[outEI, ts] = lsim(G,yis_EI,ts);
figure;
plot(ts,out,'color','m','LineWidth',1.2);hold on;
plot(ts,outZV,'-.','color','r','LineWidth',1.2);hold on;
plot(ts,outZVD,'-.','color','g','LineWidth',1.2);hold on;
plot(ts,outEI,'-.','color','b','LineWidth',1.2);hold on;
grid on;
legend('阶跃响应','ZV整形曲线响应','ZVD整形曲线响应','EI整形曲线响应');
title('欠阻尼二阶系统的振动对比')
xlabel('Time(s)');ylabel('Amp')

function y = fcnZVD(u1,u2,u3)
    b = 0.1; 
    K = exp(-b*pi/(sqrt(1-b*b))); 
    
    A1 = 1/(1+K)^2;
    A2 = 2*K/(1+K)^2;
    A3 = K^2/(1+K)^2;
    y = A1*u1 + A2*u2 + A3*u3;
end

function y = fcnEI(u1,u2,u3)
    Vexp = 0.01; %百分比
    A1 = (1+Vexp)/4;
    A2 = (1-Vexp)/2;
    A3 = (1+Vexp)/4;
    y = A1*u1 + A2*u2 + A3*u3;
end

function y = fcnZV(u1,u2)
    b = 0.1; 
    K = exp(-b*pi/(sqrt(1-b*b))); 
    
    A1 = 1/(1+K);
    A2 = K/(1+K);
    y = A1*u1 + A2*u2;
end