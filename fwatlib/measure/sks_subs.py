import numpy as np 

def splitting_intensity(Rsyn,Tsyn,dt_syn,cal_adj_source=False,si_obs = 0.,weight = 1.):
    from utils import dif1
    dRsyn = dif1(Rsyn,dt_syn)
    norm_syn = np.sum(dRsyn**2)
    norm_syn = 1. / (norm_syn + 1.0e-30)
    si_syn = -2. * np.sum(dRsyn * Tsyn) * norm_syn 

    if cal_adj_source:
        si_diff = (si_syn - si_obs) * weight
        dTcomp = dif1(Tsyn,dt_syn)
        ddRsyn = dif1(dRsyn,dt_syn)
        adjsrc_T = -2. * si_diff * dRsyn * norm_syn
        adjsrc_R = -2. * si_diff * (  
            2. * np.sum(dRsyn * Tsyn) * norm_syn**2 * ddRsyn - 
            dTcomp * norm_syn
        )
        return si_syn,adjsrc_R,adjsrc_T
    else:
        return si_syn
    