%% Creating the stream of diracs

K = 2;      % Number of Diracs
tau = 2;    % Period

[xt, loc, amp] = diracs(tau, K); % Stream of K diracs with period tau

figure;
stem(0:tau-1, xt, "filled", "MarkerSize", 5, "LineWidth", 1);	
grid on;
xlim([-0.1 tau-0.9]);
xlabel("$t$", "interpreter", "latex");
ylabel("$x(t)$", "interpreter", "latex");
title(["\bf{Original Signal}", "\bf{Stream of Diracs}", "$K=$ "+K+", $\tau=$ "+tau], "interpreter", "latex");

%% Sampling Kernel

rho = 2*K/tau; % Rate of innovation
B = rho;       % Band-width of the sinc sampling kernel ( > rho)

t = -10:0.001:10; % Continuous time
hBt = B*sinc(B*t);  % Sinc sampling kernel

figure;
plot(t, hBt, "LineWidth", 1);
grid on;
xlim([-11 11]);
xlabel("$t$", "interpreter", "latex");
ylabel("$h_B(t)$", "interpreter", "latex");
title(["\bf{Sinc Sampling Kernel}", "$h_B(t)=B\textmd{sinc}(Bt)$"], "interpreter", "latex");

%% Finding Frequency Spectrum

M = floor(B*tau/2);
N = 2*M+1;          % Number of samples > rate of innovation

% Fourier series of x(t)

Xm = zeros(1,2*M+1);
shift = M+1;
for m = -M:M
	for k = 0:K-1
		Xm(m+shift) = Xm(m+shift) + amp(k+1).*exp(-1j*2*pi*m*loc(k+1)/tau);
	end
end
Xm = (1/tau).*Xm;

% Finding the sampled yn

T = 30;
yn = zeros(1,N);
for n = 0:N-1
	shift = M+1;
	for m = -M:M
		yn(n+1) = yn(n+1) + Xm(m+shift).*exp(1j*(2*pi*m*n*T/tau));
	end
end

figure;
stem(0:N-1, real(yn), "filled", "MarkerSize", 4, "LineWidth", 1);
grid on;
xlim([0 N-1]);
xlabel("$n$", "interpreter", "latex");
ylabel("$y[n]$", "interpreter", "latex");
title("\bf{Sampled} $y_n$", "interpreter", "latex");

%% Finding time instances

% Solving the Yule-Walker System
% to find the annihilating filter coefficients

RHSX = (-1.*Xm(1+shift:K+shift))';
LHSX = zeros(K, K);
for ii = 0:K-1
	LHSX(ii+1,:) = Xm(ii+shift:-1:-K+1+ii+shift);
end

Am = [1;LHSX\RHSX];

% Plotting the zeros of annihilating filter

figure;
zplane(Am',1);

% Finding the zeros

vz = roots(Am);

% Finding the time instances of the diracs

tk = sort((tau/(2*pi)).*abs(angle(vz)));

disp("t_k");
disp(tk);

%% Finding the amplitudes of the diracs

% Solving the Vandermonde System
% to find the amplitudes of the diracs

RHSU = zeros(K, K);
for ii = 0:K-1
	RHSU(ii+1,:) = vz.^ii;
end

LHSU = Xm(0+shift:K-1+shift)';

% Finding the weights

ck = tau.*(RHSU\LHSU);

disp("c_k");
disp(ck);

%% Plotting the Result

output = zeros(1, tau);

output(tk+1)=real(ck);

figure;
stem(0:tau-1, xt, "filled", "MarkerSize", 5, "LineWidth", 1);	
grid on;
xlim([-0.1 tau-0.9]);
xlabel("$t$", "interpreter", "latex");
ylabel("$x(t)$", "interpreter", "latex");
title("\bf{Reconstruction}", "interpreter", "latex");