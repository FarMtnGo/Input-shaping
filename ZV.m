clc;clear;close all;
% 此处输入信号给定的采样时间Ts，就是输入信号的规划间隔
Ts = 1/4000;
[y1,ts] = gensig('sq',0.2,1,Ts); % 方波给定
% [u, t] = gensig(type, tau, Tf, Ts) 生成具有采样时间 Ts 的信号。
% t 从 0 到 Tf 以 Ts 的增量运行。要生成用于模拟离散时间模型的信号
% u为信号序列，t为时间序列
% type为类型，包括：sin(正弦波)，square(方波)，pulse(周期脉冲)
% tau为type类型的周期

% 方波周期是0.2s, 即5Hz, 是输入信号
% 而此处的Fs是系统自身固有频率, 是给定的系统频率，输入信息输入该系统，产生输出
Fs = 30; %Hz
wn = 2*pi*Fs; % w=2*pi*f, 系统自然频率
b = 0.1;
G = tf(wn^2,[1 2*b*wn wn^2]);
wd = wn*sqrt(1-b*b); %阻尼振荡频率
T = pi/wd; % 峰值时间

f = 1/Ts; % 采样频率,也即规划信息更新频率
tN = fix(f*T);% 取整函数

y2 = [zeros(tN,1);y1(1:end-tN,1)];
yis = [];
for n = 1:length(y1)
    yis(n,1) = fcnzv(y1(n),y2(n));
end

figure;
stairs(ts,y1,'color',[0.8500 0.3250 0.0980],'LineWidth',1.2);
hold on;
stairs(ts,yis,'-.','color',[0.4660 0.6740 0.1880],'LineWidth',1.2);
legend('方波','整形曲线')
xlabel('Time(s)');ylabel('Amp')
ylim([0 1.2])
grid on;


[y1, ts] = lsim(G,y1,ts);
[y2, ts] = lsim(G,yis,ts);
freq_domain_analysis(ts,y1)
freq_domain_analysis(ts,y2)
figure;
plot(ts,y1,'color',[0.8500 0.3250 0.0980],'LineWidth',1.2);
hold on; grid on;
plot(ts,y2,'-.','color',[0.4660 0.6740 0.1880],'LineWidth',1.2);
legend('阶跃响应','整形曲线响应');
title('欠阻尼二阶系统的振动对比')
xlabel('Time(s)');ylabel('Amp')

function y = fcnzv(u1,u2)
    b = 0.1; 
    K = exp(-b*pi/(sqrt(1-b*b))); 
    
    A1 = 1/(1+K);
    A2 = K/(1+K);
    y = A1*u1 + A2*u2;
end

