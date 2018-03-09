%pylab


xhlc = tmp1[:,1]
T1hlc = tmp1[:,4]
T2hlc = tmp2[:,4]
T3hlc = tmp3[:,4]
I1hlc = tmp1[:,2]
I2hlc = tmp2[:,2]
I3hlc = tmp3[:,2]

plot(xhlc,I1hlc,'r-')
plot(xhlc,I2hlc,'b--')
plot(xhlc,I3hlc,'g:')

xlabel('r (arcsec)', fontsize=16)
#ylabel('Core Throughput', fontsize=16)
ylabel('Mean Intensity', fontsize=16)
yscale('log')
#ylim(0.005,0.040)
ylim(1e-12,1e-10)
title('HLC')
tick_params(axis='both', which='major', labelsize=14)
savefig('HLC_I', dpi=300, transparent=True)

