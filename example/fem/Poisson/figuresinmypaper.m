%% figures in mypaper
%% 求解区域
figure;
slope = 0.1;

x = 0:0.01:1;
ybot = -0.1*x;
ytop = ybot + 1;

subplot(2,3,1);
hold on;
plot(x,ybot);
plot(x,ytop);
