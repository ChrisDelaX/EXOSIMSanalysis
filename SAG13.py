import numpy as np
import astropy.units as u

# numerical integration of parametric occurrence rate distribution

def calc_maxintTime(self, TL):


# bin definition
Rrange = [4200, 6400]*u.km
Prange = [320, 640]*u.day

# reshape input
Rrange = np.array(Rrange.to('earthRad'),ndmin=1)
Prange = np.array(Prange.to('year'),ndmin=1)


# distribution parameters



# numerical sampling grid
dlogR = 0.001
dlogP = 0.001
logR = log(Rmin):dlogR:log(Rmax); # in Earth radii
logP = log(Pmin):dlogP:log(Pmax); # in years
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


###

def occur_density(Rp, period):
    
    Gamma = [0.38, 0.73]
    alpha  = [-0.19, -1.18]
    beta   = [0.26, 0.59]
    Rbreak = 3.4*u.earthRad
    
    R = Rp.to('earthRad').value
    P = period.to('year').value
    dN = np.zeros(len(R))
    mask = R < Rbreak
    dN[mask]  = Gamma[0] * R**alpha[0] * P**beta[0]
    dN[~mask] = Gamma[1] * R**alpha[1] * P**beta[1]
    
    return -np.dot(diff,np.dot(icov,diff))/2.0













