import os.path,EXOSIMS,EXOSIMS.MissionSim
scriptfile = os.path.join(EXOSIMS.__path__[0],'Scripts','script_OptSys.json')
sim = EXOSIMS.MissionSim.MissionSim(scriptfile)
%pylab

#nplan = 1e6
#s, dmag = sim.Completeness.genplans(nplan)
#plot(s,dmag, 'bo')     # Testing completeness

###

import os.path,json,EXOSIMS.Prototypes.OpticalSystem
scriptfile = os.path.join(EXOSIMS.__path__[0],'Scripts','script_OptSys.json')
specs = json.loads(open(scriptfile).read())


import numpy as np
OS = sim.OpticalSystem

a = 1
if a==1:
    print 1
elif a==2:
    print 2
else:
    print 3

dict = {'feh':1}
a = dict.get('feh',2)
a = dict.get('notfeh',raise ValueError("No imager defined."))


# solar longitude of 135 degrees, where the local zodi is minimized, minor solar term "Autumn Commences"
lwz = 3;
fsz = 24;
grid('on');xlabel('ecliptic latitude sin($\beta$)',fontsize=1.5*fsz);ylabel('Variation of zodiacal light $f_{\beta}$',fontsize=1.5*fsz)
plot(sinbeta,f2,label='Leinert table')
plot(sinbeta,f1,label='Stark empiric')
plot(sinbeta,f3,label='Lindler empiric')
legend(fancybox=True, framealpha=0.3, loc='lower right', fontsize=fsz);

sinbeta=np.linspace(0,np.sin(np.radians(ss_elongation)),21)

beta = np.linspace(25,50,6)
sinbeta = np.sin(np.radians(beta))

# Αα    Alpha     Νν    Nu
# Ββ    Beta      Ξξ    Xi
# Γγ    Gamma     Οο    Omicron
# Δδ    Delta     Ππ    Pi
# Εε    Epsilon   Ρρ    Rho
# Ζζ    Zeta      Σσς   Sigma
# Ηη    Eta       Ττ    Tau
# Θθ    Theta     Υυ    Upsilon
# Ιι    Iota      Φφ    Phi
# Κκ    Kappa     Χχ    Chi
# Λλ    Lambda    Ψψ    Psi
# Μμ    Mu        Ωω    Omega

