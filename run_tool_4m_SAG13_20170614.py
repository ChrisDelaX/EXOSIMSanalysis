import EXOSIMS, EXOSIMS.MissionSim
#from EXOSIMS.util.run_one import run_one
import os.path,json
import datetime


# JSON script file
filename = 'HabEx_4m_SAG13_cluster_20170613.json'
# Path to JSON script file
scriptpath = '/proj/exep/rhonda/HabEx/SAG13'
#scriptpath = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/testing/Scripts')
scriptfolder = os.path.normpath(os.path.expandvars(scriptpath))
# Path to saved results
savepath = '/proj/exep/rhonda/HabEx/SAG13/out'
#savepath = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/testing/results')
savepath = os.path.join(savepath, filename[:-5]) #removes '.json'
savefolder = os.path.normpath(os.path.expandvars(savepath))
# Check if save folder exists and create if necessary
if not os.path.exists(savefolder):
    os.makedirs(savefolder)

# Number of simulations to run
nbsim = 2 #must be an integer
# Selected actions:
DO_ipyparallel = True
DO_run_ensemble = False
DO_load_results = False

###############################

# DO_run_ensemble: run an ensemble of nb_sim simulations, and store results in save folder
if DO_run_ensemble is True:
    scriptfile = os.path.join(scriptfolder,filename)
    for i in xrange(int(nbsim)):
        print '\nSimulation number ', i+1, '/', int(nbsim)
        reload(EXOSIMS)
        reload(EXOSIMS.MissionSim)
        sim = EXOSIMS.MissionSim.MissionSim(scriptfile)
        SS = sim.SurveySimulation
        SS.run_sim()
        SS._outspec['DRM'] = SS.DRM
        spcpath = os.path.join(savefolder, str(sim._outspec['seed']))
        sim.genOutSpec(spcpath)

if DO_ipyparallel is True:
    # define the run_one function (you can specify 'genNewPlanets' and 'rewindPlanets')
    def run_one(): 
        # wrap the run_sim in a try/except loop
        nbmax = 10
        for attempt in range(nbmax):
            try:
                # run one survey simulation
                SS.run_sim()
                res = SS.DRM[:]
            except:
                # if anything goes wrong, reset simulation
                sim.reset_sim()
            else:
                break
        else:
            raise ValueError("Unsuccessful run_sim after %s reset_sim attempts"%nbmax)
        
        # reset simulation at the end of each simulation
        SS.reset_sim(genNewPlanets=True, rewindPlanets=True)
        
        # need to pass savefolder to run_one
        import os.path
        filename = 'HabEx_4m_SAG13_cluster_20170613.json'
        savepath = '/proj/exep/rhonda/HabEx/SAG13/out'
#        savepath = os.path.expandvars('$HOME/CODE/Python/EXOSIMS/testing/results')
        savepath = os.path.join(savepath, filename[:-5]) #removes '.json'
        savefolder = os.path.normpath(os.path.expandvars(savepath))
        if not os.path.exists(savefolder):
           os.makedirs(savefolder)
        # store DRM of each survey
        import pickle
        from random import randint
        pklname = ''.join(["%s" % randint(0, 9) for num in range(9)]) + '.pkl'
        pklpath = os.path.join(savefolder, pklname)
        with open(pklpath, 'wb') as f:
            pickle.dump(res, f)
        # store OutSpec of each survey
        spcname = pklname[:-4] + '.spc'
        spcpath = os.path.join(savefolder, spcname)
        SS.genOutSpec(spcpath)
        
        return res
    
    scriptfile = os.path.join(scriptfolder, filename)
    sim = EXOSIMS.MissionSim.MissionSim(scriptfile)
    ens = sim.run_ensemble(nbsim, run_one=run_one)
    timenow = datetime.datetime.now()
    format = "%Y%m%d_%H%M"
#    sim.genOutSpec(savefolder + '/' + timenow.strftime(format))


# DO_load_results: load the results stored in save folder and build 
if DO_load_results is True:
    names = [x for x in os.listdir(savefolder) if x not in '.DS_Store']
    nb_load_sim = len(names)
    ENS_obs = []
    ENS_det = []
    ENS_md = []
    ENS_outside_iwa = []
    ENS_outside_owa = []
    ENS_char = []
    ENS_uniq = []
    for name in names:
        specs = json.loads(open(savefolder+'/'+name).read())
        DRM = specs['DRM']
        nobs = len(DRM)
        ENS_obs.append(nobs)
        det,md,outside_iwa,outside_owa = [],[],[],[]
        for i in xrange(nobs):
            det.append(len([x for x in DRM[i]['det_status'] if x == 1]))
            md.append(len([x for x in DRM[i]['det_status'] if x == 0]))
            outside_iwa.append(len([x for x in DRM[i]['det_status'] if x == -1]))
            outside_owa.append(len([x for x in DRM[i]['det_status'] if x == -2]))
        ENS_det.append(sum(det))
        ENS_md.append(sum(md))
        ENS_outside_iwa.append(sum(outside_iwa))
        ENS_outside_owa.append(sum(outside_owa))
