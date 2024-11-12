clc;clear;close all;
% s曲线仿真
S_Trajdata = readtable('PosRef_Traj.xlsx');
S_Trajdata = table2array(S_Trajdata);
Fhz = 4000; % 采样频率
NN = size(S_Trajdata,1);
Ts = 0:1/Fhz:(NN-1)/Fhz;

% ZV 两脉冲 脉冲时滞 t1 = 0; t2 = T;
wn = 32; 
b = 0.2;
G1 = tf(wn*wn,[1 2*wn*b wn*wn]); % 原系统
step(G1)

K = exp(-b*pi/(sqrt(1-b*b))); 
T = pi/((wn*sqrt(1-b*b))); % 峰值时间
A1 = 1/(1+K);
A2 = K/(1+K);
um = S_Trajdata(:,2); % 轨迹
t = fix(T*Fhz);
% 输入整形
% u1(:,1) = A1*um; % 第一个脉冲
% u2(:,1) = A2*um; % 第二个脉冲
% u1(t+1:NN,1) = u1(t+1:NN,1)+ u2(1:NN-t,1); % 延迟
% 另一种写法
% u1(1:t,1) = A1*um(1:t,1);
% u1(t+1:NN,1) = A1*um(t+1:NN,1) + A2*um(1:NN-t,1);
% 另一种写法
um_ = [zeros(t,1);um(1:NN-t,1)];
u1 = A1*um + A2*um_;

figure;
y1 = lsim(G1,u1,Ts);
plot(Ts,y1);
title('1')
figure;
plot(Ts,um,Ts,u1);
legend('未加输入整形','输入整形');