class FWAT_MODEL_CONST:
    def __init__(self,mdtype='iso',use_model_set=1) -> None:
        # iso model by default, we use parameter set vp_vs_rho
        # = 0 kappa_mu_rho
        # = 1  vp_vs_rho
        # = 2 vpvs_vs_rho  
        self.MODEL_TYPE = mdtype
        self.USE_MODEL_SET = use_model_set

        if mdtype != 'iso':
            print(f"not implemented for modeltype = {mdtype}")
            exit(1)

def convert_md(model,backward=False):
    '''
    convert base parameters to another model set
    # base parameters are vp/vs/rho for iso model
    '''
    # make a copy 
    md1 = model.copy()

    # constants
    p = FWAT_MODEL_CONST()

    if p.MODEL_TYPE == 'iso':
        # base parameters are vp/vs/rho
        if p.USE_MODEL_SET == 2:
            if backward:
                md1[0,...] = model[0,...] * model[1,...] 
            else:
                md1[0,...] = model[0,...] / model[1,...]  
    
    return md1 

def convert_kernel(md,md_kl):
    """
    convert kernel
    # base parameters are vp/vs/rho for iso model
    """
    md1_kl = md_kl.copy()

    # constants
    p = FWAT_MODEL_CONST()

    if p.MODEL_TYPE == 'iso':

        if p.USE_MODEL_SET == 2: # vpvs_vs_rho
            md1_kl[0,...] = md_kl[0,...].copy()
            md1_kl[1,...] = md_kl[0,...] + md_kl[1,...] 
        elif p.USE_MODEL_SET == 0: # kappa_mu_rho
            # please change it 
            md1_kl = md * md_kl

        return md1_kl

def get_gradname_list(smooth=False):
    p = FWAT_MODEL_CONST()

    if p.MODEL_TYPE == 'iso':
        grad_list = ['alpha_kernel','beta_kernel','rhop_kernel']
    
    if smooth:
        for g in grad_list:
            g = g + '_smooth'
    
    return grad_list

def get_model_name_list():
    p = FWAT_MODEL_CONST()
    if p.MODEL_TYPE == 'iso':
        mod_list = ['vp','vs','rho']
    
    return mod_list

def get_direc_name_list():
    p = FWAT_MODEL_CONST()
    if p.MODEL_TYPE == 'iso':
        direc_list = ['dbulk','dbeta','drho']
        #grad_list = ['alpha_kernel_smooth','beta_kernel_smooth','rhop_kernel_smooth']
        if p.USE_MODEL_SET == 2:
            direc_list[0] = 'dvpvs'
    
    return direc_list
