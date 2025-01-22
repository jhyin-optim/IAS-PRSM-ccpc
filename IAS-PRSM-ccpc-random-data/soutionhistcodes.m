figure
subplot(3,3,1)
plot(y1,'b','LineWidth',1.5)
ylim([-20 20]), 
ylabel('x^*'); 
legend('AS-ADMM')
grid on

ylabel('x^*'); 
subplot(3,3,2)
plot(y2,'r','LineWidth',1.5)
ylim([-20 20]), 
ylabel('x^*'); 
legend('SAS-ADMM')
grid on

subplot(3,3,3)
plot(y3,'g','LineWidth',1.5)
ylim([-20 20]), 
ylabel('x^*'); 
legend('IASPRSM-ccpc')
grid on
subplot(3,3,[4 7])
[hc,h] = hist(y1,[-.1:.01:.1]);
bar(h,hc,'b')
axis([-.1 .1 -10 800])
ylabel('hist(x^*)'); 
legend('AS-ADMM')
grid on

subplot(3,3,[5 8])
[hc,h] = hist(y2,[-.1:.01:.1]);
bar(h,hc,'r')
axis([-.1 .1 -10 800])
ylabel('hist(x^*)'); 
legend('SAS-ADMM')
grid on

subplot(3,3,[6 9])
[hc,h] = hist(y2,[-.1:.01:.1]);
bar(h,hc,'g')
axis([-.1 .1 -10 800])
ylabel('hist(x^*)'); 
legend('IASPRSM-ccpc')
grid on
% Comparison of the solutions obtained by AS-ADMM (blue, left) and by las-ADMM 
% (red, right) for Problem (5.1) on the mnist dataset (120s).