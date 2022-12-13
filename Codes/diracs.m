function [signal, loc, amp] = diracs(tau,K)

signal = zeros(1,tau);
% loc = sort(randperm(tau,K));
loc = [0, 1];
% amp = 2*rand(1,K)-1;
amp = [1, -2];
signal(loc+1) = amp;

end