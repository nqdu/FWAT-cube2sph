import numpy as np
from fwat.optimize.model import FwatModel



def user_defined_misfit(weight,cmodel:np.ndarray):
    """
    Example of user-defined misfit function.
    This function computes a simple L2 norm of the model parameters
    as a demonstration. Users can modify this function to implement
    their own misfit calculations based on specific requirements.

    Args:
        cmodel: np.ndarray
            Array of model parameters.

    Returns:
        float
            Computed misfit value.
    """
    # Simple L2 norm as an example misfit
    misfit = np.sum(weight * cmodel**2) * 0.5
    kl = weight * cmodel
    return misfit,kl

def _Index(m,n):
    assert(m<=n)
    idx = m * 6 + n - (m * (m + 1)) // 2

    return idx

def get_base_model(mdtype,kltype):
    if mdtype  =='dtti':
        vsh = 3.5
        vsv = 3.5

        if kltype > 1:
            vsh = 3.7
        
        vph = 1.732 * vsv
        vpv = 1.732 * vsh
        rho = 2.8 
        eta = 1.0
        A = rho * vph**2 
        L = rho * vsv**2 
        C = rho * vpv**2 
        N = rho * vsh**2
        F = eta * (A - 2* L)

        G0 = 0.05 
        phi = np.pi / 3 
        gcp = G0 * np.cos(2 * phi)
        gsp = G0 * np.sin(2 * phi)
        gc = gcp * L
        gs = gsp * L

        cmodel = np.zeros((22,1),dtype='f4')
        cmodel[_Index(0,0),...] = A
        cmodel[_Index(0,1),...] = A - 2 * N  
        cmodel[_Index(0,2),...] = F
        cmodel[_Index(1,1),...] = A
        cmodel[_Index(1,2),...] = F
        cmodel[_Index(2,2),...] = C
        cmodel[_Index(3,3),...] = L - gc
        cmodel[_Index(3,4),...] = -gs
        cmodel[_Index(4,4),...] = L + gc
        cmodel[_Index(5,5),...] = N
        cmodel[-1,...] = rho
    else:
        print("not implemented mdtype")
        exit(1)

    return cmodel

def main():
    # kltype and model type test
    test_cases_dict = {
        "dtti":[1,2,3]
    }
    # fix random seed
    np.random.seed(42)
    
    for mtype, kltypes in test_cases_dict.items():
        for kltype in kltypes:
            
            # generate model 
            cmodel = get_base_model(mtype,kltype)
            weight = np.random.rand(cmodel.shape[0],1)

            # compute misfit and kernels for base model
            misfit,kl_base = user_defined_misfit(weight,cmodel)

            M = FwatModel(filename=None,mdtype = mtype,kltype=kltype)
            
            # convert from base model to defined model
            model = M.convert_model(cmodel)
            _,kl_analytic = M.convert_kl(cmodel,kl_base)

            # use FD to compute kernel
            kl_fd = np.zeros_like(kl_analytic)
            delta = 1.0e-2
            for i in range(model.shape[0]):
                model_p = model.copy()
                model_m = model.copy()
                model_p[i,0] *= (1. + delta)
                model_m[i,0] *= (1. - delta)

                # compute misfit
                cmodel_p = M.convert_model(model_p,backward=True)
                cmodel_m = M.convert_model(model_m,backward=True)
                misfit_p,_ = user_defined_misfit(weight,cmodel_p)
                misfit_m,_ = user_defined_misfit(weight,cmodel_m)

                # FD kernel
                kl_fd[i,0] = (misfit_p - misfit_m) / (2. * delta * model[i,0]) 
            
            # change to relative change
            if mtype == 'dtti' and kltype == 3:
                kl_fd[:3,...] *= model[:3,...]
            if mtype == 'dtti' and kltype == 1:
                kl_fd[:3,...] *= model[:3,...]
            
            if mtype == 'dtti' and kltype == 2:
                kl_fd[0:5,...] *= model[0:5,...]
            
            # check
            diff = np.abs(kl_fd - kl_analytic) / (np.abs(kl_analytic) + 1.e-20) * 100
            direc_names = M.get_direc_names()
            # rel_err = np.max(diff) / (np.max(np.abs(kl_analytic)) + 1.e-20)
            # print(f"Model type: {mtype}, KL type: {kltype}, Relative error: {rel_err:.3e}")
            # print(kl_fd - kl_analytic)
            for i in range(kl_analytic.shape[0]):
                print(f"Model type: {mtype}, KL type: {kltype}, Direction: {direc_names[i]}, kl_fd: {kl_fd[i,0]:.6e}, kl_analytic: {kl_analytic[i,0]:.6e}, relative error: {diff[i,0]}%")
            print("")

if __name__ == "__main__":
    main()