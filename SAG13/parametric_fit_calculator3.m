clear all;
% numerical integration of parametric occurrence rate distribution

% distribution parameters
Gamma  = [0.38 0.73];
alpha  = [-0.19 -1.18];
beta   = [0.26 0.59];
Rbreak = [0 3.4 Inf];

% bin definition
Rmin = .667;
Rmax = 3.4;
Pmin = 320/365;
Pmax = 640/365;

% numerical sampling grid
dlogR = 0.001;
dlogP = 0.001;
logR = log(Rmin):dlogR:log(Rmax); % in Earth size
logP = log(Pmin):dlogP:log(Pmax); % in years
R = exp(logR);
P = exp(logP);
[logPP logRR] = meshgrid(logP, logR);
[PP RR] = meshgrid(exp(logP), exp(logR));

% distribution across sampling grid 
n = zeros(size(RR));
for i = 1:length(Gamma)
    region = (Rbreak(i) <= RR).*(RR < Rbreak(i+1));
    n = n + Gamma(i) * RR.^alpha(i) .* PP.^beta(i) .* region;
end

% Integrated occurrence in the bin
N = trapz(trapz(n.* dlogR * dlogP))