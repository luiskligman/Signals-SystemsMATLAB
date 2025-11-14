clear all;
close all;
clc;

%% Solving for a0, an, bn and stating xN
T0 = 2;  % Period
w0 = pi;

% x(t) = -1 from 0 to 1
%         1 from 1 to 2

% Calculated a0, an, and bn on paper due to the complexity of solving them
% in matlab
a0 = 0;
an = 0;
bN = @(N) (-4./(N*pi)).*(mod(N,2)==1);  % odd n -> -4/npi, even n -> 0

t = linspace(-5,5,1000);

%% Question 1: cosine/sin Fourier Series Representation
xN = @(t,N) sum(bN(1:N).' .* sin((1:N)'*w0.*t), 1);

for N = [1, 3, 5, 10, 100]
figure(1)
plot(t,xN(t,N),'DisplayName',['xN_{' num2str(N) '}(t)'])
hold on
grid on
xlabel('time')
ylabel('xN(t)')
end
legend('show')
title('Cosine / Sin Fourier Series Representation')


%% Question 2: amplitude/phase Fourier Series Representation
% c_n = sqrt( (b_n)^2 )
% theta_n = -tan^-1 (b_n/0) = -tan^-1 (-inf) --> pi/2

C_n = @(N) abs((-4./(N*pi)) .* (mod(N,2)==1));
theta_n = pi/2;

xN_amp = @(t,N) sum( C_n(1:N).' .* cos((1:N)'*w0.*t + theta_n), 1 );

for N = [1, 3, 5, 10, 100]
figure(2)
plot(t,xN_amp(t,N),'DisplayName',['xN amp_{' num2str(N) '}(t)'])
hold on
grid on
xlabel('time')
ylabel('xN amp(t)')
end
legend('show')
title('Amplitude / Phase Fourier Series Representation')


%% Question 3: complex-exponential Fourier Series Representation

% xn = c0 when n = 0
% xn = cn/2 * e^jtheta_n when n > 0
%x_n = x_n^*

c_of = @(n) (mod(abs(n),2)==1).* (2j./(n*pi));  % 0 for even or n=0, 2j/(n*p) for od
xN_exp = @(t,N) real( sum( c_of(-N:N).' .* exp(1j* (-N:N)' *w0 .* t), 1));


hold on;
grid on;
for N = [1, 3, 5, 10, 100]
figure(3)
plot(t,xN_exp(t,N),'DisplayName',['xN exp_{' num2str(N) '}(t)'])
hold on
grid on
xlabel('time')
ylabel('xN exp(t)')
end
legend('show')
title('Complex-Exponential Fourier Series Representation')
