[y,ts] = gensig('sq',0.2,1,1/4000); % 方波给定
yis = [];

for n = 1:length(y)
    yis(n,1) = fcnzv(y(n),40);
end

figure;
stairs(ts,y,'color',[0.8500 0.3250 0.0980],'LineWidth',1.2);
hold on;
stairs(ts,yis,'-.','color',[0.4660 0.6740 0.1880],'LineWidth',1.2);
legend('方波','整形曲线')
xlabel('Time(s)');ylabel('Amp')
ylim([0 1.2])
grid on;

wn = 40*2*pi;
b = 0.1;
G = tf(wn^2,[1 2*b*wn wn^2]);
[y1, ts] = lsim(G,y,ts);
[y2, ts] = lsim(G,yis,ts);
figure;
plot(ts,y1,'color',[0.8500 0.3250 0.0980],'LineWidth',1.2);
hold on; grid on;
plot(ts,y2,'-.','color',[0.4660 0.6740 0.1880],'LineWidth',1.2);
legend('阶跃响应','整形曲线响应');
title('欠阻尼二阶系统的振动对比')
xlabel('Time(s)');ylabel('Amp')

function y = fcnzv(u,fs)
    persistent u_k_1   array_u 
    b = 0.1; wn = fs*2*pi;
    K = exp(-b*pi/(sqrt(1-b*b))); 
    T = pi/((wn*sqrt(1-b*b))); % 峰值时间
    f = 4000; % 采样频率
    t2 = fix(f*T);
    A1 = 1/(1+K);
    A2 = K/(1+K);
    if isempty(u_k_1); u_k_1 = 0;array_u = zeros(2000,1);end
    y = A1*u + A2*array_u(t2);
    for i = 2000:-1:2
        array_u(i) = array_u(i-1);
    end
    array_u(1) = u;
end