% 输入整形 仿真
wn = 30; % 无阻尼固有频率
b = 0.2; % 阻尼比
G1 = tf(wn*wn,[1 2*wn*b wn*wn]); % 原系统
step(G1)

wn = 32; % 模型失配参数
b = 0.2;
K = exp(-b*pi/(sqrt(1-b*b))); 
T = pi/((wn*sqrt(1-b*b))); % 峰值时间

Index = 15; % 脉冲时滞 时间因子
TN = 500; % 一个脉冲时滞 T 线条点数
TNm = Index*TN; % 线条总点数
Ts = linspace(0,Index*T,TNm)';

% ZV 两脉冲 脉冲时滞 t1 = 0; t2 = T;
A1 = 1/(1+K);
A2 = K/(1+K);

u1(1:TN,1) = A1; % 第一个脉冲
u1(TN+1:TNm,1) = A2+A1; % 第二个脉冲
y1 = lsim(G1,u1,Ts);
y2 = step(G1,Ts);
plot(Ts,u1,Ts,y2,Ts,y1);
title('阶跃响应');legend('ZV输入整形脉冲','未加输入整形','输入整形');
axis([0 1 0 1.6]);

% EI 脉冲时滞 t1 = 0; t2 = T; t3 = 2T;
Vexp = 0.1;
Ae1 = (1+Vexp)/4;
Ae2 = (1-Vexp)/2;
Ae3 = (1+Vexp)/4;
u2(1:TN,1) = Ae1;
u2(TN+1:2*TN,1) = Ae1 +Ae2;
u2(2*TN+1:TNm,1) = Ae1 +Ae2 + Ae3;
ye1 = lsim(G1,u2,Ts);
ye2 = step(G1,Ts);
plot(Ts,u2,Ts,ye2,Ts,ye1);
title('阶跃响应');legend('EI输入整形脉冲','未加输入整形','输入整形');
axis([0 1 0 1.6]);
