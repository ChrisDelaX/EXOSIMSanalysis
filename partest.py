import numpy as np
import EXOSIMS
import EXOSIMS.MissionSim

def run_one():
    SS.run_sim()
    res = SS.DRM[:]
    SS.reset_sim()
    return res


if __name__ == "__main__":
    sim = EXOSIMS.MissionSim.MissionSim('parscript.json')
    res = sim.run_ensemble(8,run_one=run_one)

    print res

