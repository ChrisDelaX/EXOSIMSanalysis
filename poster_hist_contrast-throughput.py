# fsz = 34
# figure;grid('on');
# xlabel('max separation (arcsec)',fontsize=fsz);ylabel('Number of star targets', fontsize=fsz);
# title('Before filtering',fontsize=fsz)
# ymax = 1200
# hist(TL.WA_before,bins=50);xlim(0,0.45);ylim(0,ymax)
# plot([OS.IWA.value,OS.IWA.value],[0,ymax],'k--',linewidth=lwz)
# xticks([0.15,0.30,0.45], fontsize = fsz)
# yticks([0,400,800,ymax], fontsize = fsz)
# tight_layout()
# savefig('before.png', dpi=300, transparent=True);
# 
# figure;grid('on');
# xlabel('max separation (arcsec)',fontsize=fsz);ylabel('Number of star targets', fontsize=fsz);
# title('After filtering',fontsize=fsz)
# ymax = 60
# hist(TL.WA_after,bins=50);xlim(0,0.45);ylim(0,ymax)
# plot([OS.IWA.value,OS.IWA.value],[0,ymax],'k--',linewidth=lwz)
# xticks([0.15,0.30,0.45], fontsize = fsz)
# yticks([0,20,40,ymax], fontsize = fsz)
# tight_layout()
# savefig('after.png', dpi=300, transparent=True);



# lwz = 3;
# fsz = 25
# import astropy.io.fits as fits
# import scipy.interpolate
# import scipy.optimize
# Tpth = os.path.normpath(os.path.expandvars('$HOME/INSTRUMENTS/AFTA-CGI/hlc20140623-139/hlc_20140623-139_polx_0.4mas_jitter_10x_throughput.fits'))
# with fits.open(Tpth) as tmp:
#     dat = tmp[0].data
#     WAT = dat[0]
#     T = dat[1]
# Tinterp = scipy.interpolate.interp1d(WAT, T, kind='cubic',fill_value=np.nan, bounds_error=False)
# Tmax = scipy.optimize.minimize(lambda x:-Tinterp(x),WAT[np.argmax(T)],bounds=((np.min(WAT),np.max(WAT)),) )
# Tmax = -Tmax.fun[0]
# Cpth = os.path.normpath(os.path.expandvars('$HOME/INSTRUMENTS/AFTA-CGI/hlc20140623-139/hlc_20140623-139_polx_0.4mas_jitter_10x_contrast.fits'))
# with fits.open(Cpth) as tmp:
#     dat = tmp[0].data
#     WAC = dat[0]
#     C = dat[1]
# f=figure(1);
# #f.suptitle('HLC, jitter=0.4mas, f_PP=0.1', fontsize = fsz, y=1.05)
# f.subplots_adjust(hspace=0)
# ax1 = subplot2grid((2,2), (0, 0),colspan=2)
# ax2 = subplot2grid((2,2), (1, 0),colspan=2);ax2.set_yscale('log');ax2.set_xlim(0,.5);ax2.set_ylim(1e-10,2e-8);
# ax1.plot(WAT,T,linewidth=lwz);ax1.grid('on')
# ax1.plot([IWA,IWA],[0,Tmax/2],'k--',linewidth=lwz);ax1.plot([OWA,OWA],[0,Tmax/2],'k--',linewidth=lwz);ax1.plot([0,OWA],[Tmax/2,Tmax/2],'k--',linewidth=lwz);
# ax2.plot(WAC,C,linewidth=lwz);ax2.grid('on')
# ax2.plot([IWA,IWA],[0,2e-8],'k--',linewidth=lwz);ax2.plot([OWA,OWA],[0,2e-8],'k--',linewidth=lwz);
# ax1.set_yticks([0.004,0.008,0.012]);ax1.set_xticks([0,.1,.2,.3,.4,.5]);
# ax1.set_xticklabels([])
# ax2.set_yticks([1e-10,1e-9,1e-8]);ax2.set_xticks([0,.1,.2,.3,.4,.5]);
# ax1.set_ylabel('Throughput',fontsize=1.3*fsz);ax2.set_ylabel('Contrast',fontsize=1.3*fsz);
# ax2.set_xlabel('Separation (arcsec)',fontsize=1.3*fsz)
# ax1.tick_params(labelsize=fsz)
# ax2.tick_params(labelsize=fsz)
# tight_layout()
# savefig('TandC.png', dpi=300, transparent=True);