figure
subplot(3,2,1)
plot(y1,'b','LineWidth',1.5)
ylim([-20 20]), 
ylabel('x^*'); 
legend('AS-ADMM')
grid on
ylabel('x^*'); 
subplot(3,2,2)
plot(y2,'r','LineWidth',1.5)
ylim([-20 20]), 
ylabel('x^*'); 
legend('SAS-ADMM')
grid on
subplot(3,2,[3 5])
[hc,h] = hist(y1,[-.1:.01:.1]);
bar(h,hc,'b')
axis([-.1 .1 -10 800])
ylabel('hist(x^*)'); 
legend('AS-ADMM')
grid on

subplot(3,2,[4 6])
[hc,h] = hist(y2,[-.1:.01:.1]);
bar(h,hc,'r')
axis([-.1 .1 -10 800])
ylabel('hist(x^*)'); 
legend('SAS-ADMM')
grid on
% Comparison of the solutions obtained by AS-ADMM (blue, left) and by las-ADMM 
% (red, right) for Problem (5.1) on the mnist dataset (120s).