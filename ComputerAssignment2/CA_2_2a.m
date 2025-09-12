% ELCT 222 --- CA2-2a (minimal)
clear; clc; close all;

A = 1; B = 1; T1 = 2; T2 = 1;          % set the case you need

yfun = @(t) A*B*max(0, min(T1, t) - max(0, t - T2));  % analytic result

t = linspace(-0.5, T1+T2+0.5, 800);
x = A*(t>=0 & t<T1);
h = B*(t>=0 & t<T2);
y = yfun(t);

plot(t,x,'LineWidth',1.5); hold on;
plot(t,h,'LineWidth',1.5);
plot(t,y,'LineWidth',2);
grid on; xlabel('t (s)'); ylabel('Amplitude');
legend('x(t)','h(t)','y(t)=x*h','Location','best');
title('Q2(a): Convolution of two rectangular pulses');