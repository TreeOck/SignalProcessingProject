% Continuous Time
T0 = 0.001;
tcont = -10:T0:10;

% Original Function
f = 1/5;
xt = sin(2*pi*f*tcont);

% Discrete Time
T = 0.1;
tsamples = -10:T:10;

% Sampled Function
N = length(tsamples);
xn = zeros(N, 1);
for k = 1:N
	xn(k) = xt((k-1)*(T/T0) + 1);
end

n = -10/T:10/T;
y = sinc_reconstruction(n, xn, T, tcont);

% Test
figure;
sgtitle("\bf{Sinc Reconstruction}", "interpreter", "latex");

subplot(3,1,1);
plot(tcont,xt,"b");
grid on;
title("\bf{Original Signal}", "interpreter", "latex");
xlabel("$t$", "interpreter", "latex");
ylabel("$x(t)$", "interpreter", "latex");

subplot(3,1,2);
stem(tsamples,xn,"filled","MarkerSize",2,"Color","r");
grid on;
title("\bf{Sampled Signal}", "interpreter", "latex");
xlabel("$t$", "interpreter", "latex");
ylabel("$x(nT_s)$", "interpreter", "latex");

subplot(3,1,3);
plot(tcont,y,"r","LineWidth",1.5);
hold on;
plot(tcont,xt,"b","LineWidth",2,"LineStyle","--");
hold off;
grid on;
legend("Sinc Reconstruction","Original Signal");
title("\bf{Reconstruction}", "interpreter", "latex");
xlabel("$t$", "interpreter", "latex");
ylabel("Signal", "interpreter", "latex");