function y = sinc_reconstruction(n,xn,T,tcont)

	% n     - the integer locations of the samples x[n]
	% xn    - the sampled signal x[n] = x(n*Ts)
	% T     - the sampling interval
	% tcont - the time-grid for reconstruction of xr
	% y     - the reconstructed signal over the time-grid tcont

	ws = pi/T;

	lent = length(tcont);
	lens = length(n);

	y = zeros(1, lent);

	for t = 1:lent
		% For every time instant 't'
		% We calculate the sinc-interpolation over all samples
		% And add it to the running sum
		% xr(t) = xr(t) + Ts * x(n*Ts) * (wc/pi) * sinc(wc*(t-n*Ts)/pi)
		for s = 1:lens
			y(t) = y(t) + T*xn(s)*(ws/pi)*sinc(ws*(tcont(t)-n(s)*T)/pi);
		end
	end
end