from typing import Optional
import numpy as np 
from fwat.const import PARAM_FILE
from .dtti_generated import KERNEL_CONVERTERS, MODEL_CONVERTERS

def _Index(m,n):
    assert(m<=n)
    idx = m * 6 + n - (m * (m + 1)) // 2

    return idx

class FwatModel:
    def initialize(self,mdtype='iso',kltype=1) -> None:
        """
        Initialize the FwatModel.

        Parameters
        ----------
        mdtype : str
            Model type. Supported values:
            - 'iso'  : Isotropic
            - 'dtti' : Transversely isotropic model 
        kltype : int
            Kernel type. Supported values:
            For iso:
                - 1 : vp, vs, rho
                - 2 : vp/vs, vs, rho

            For dtti:
                * 1 : vp, vs, rho, gcp, gsp
                * 2 : vph, vpv, vsh, vsv, rho, eta, gcp, gsp
                * 3 : vp,vs,rho,(vph-vpv)/vpv,(vsh-vsv)/vsv, eta, gcp, gsp
        """
        self._mdtype = mdtype
        self._kltype = kltype

        if mdtype not in ['iso','dtti']:
            raise NotImplementedError(f"Model type '{mdtype}' is not implemented. Supported types are 'iso' and 'dtti'.")

    def __init__(self,filename:Optional[str]=PARAM_FILE,mdtype='iso',kltype=1) -> None:
        """
        when filename is provided, mdtype and kltype will be overwritten
        """
        self._mask_vars = []
        if filename is not None:
            import yaml
            with open(filename,"r") as f:
                pdict = yaml.safe_load(f)['optimize']
            mdtype = pdict['MODEL_TYPE']
            kltype = pdict['KERNEL_SET']

            # mask option
            
            if 'MASK_VARS' in pdict.keys():
                self._mask_vars = pdict['MASK_VARS']
    
        self.initialize(mdtype,kltype)

    def model_names(self):
        if self._mdtype == 'iso':
            mod_list = ['vp','vs','rho']
        else:
            mod_list = [f'c{i+1}{j+1}' for i in range(6) for j in range(i,6)]
            mod_list += ['rho']
        
        return mod_list
    
    def direc_names(self):
        if self._mdtype == 'iso':
            direc_list = ['dalpha','dbeta','drho']
            if self._kltype == 2:
                direc_list[0] = 'dvpvs'
        elif self._mdtype == "dtti":
            direc_list = [f'dc{i+1}{j+1}' for i in range(6) for j in range(i,6)]
            direc_list += ['drho']
            if self._kltype == 1:
                # vp,vs,rho,gc_nodim,gs_nodim  Zhu et al 2015, GJI, (25,26) 
                direc_list = ["dalpha","dbeta","drho","dGcp","dGsp"]
            elif self._kltype == 2:
                direc_list = ["dalphah","dalphav","dbetah","dbetav","drho","deta","dGcp","dGsp"]
            elif self._kltype == 3:
                direc_list = ["dalpha","dbeta","drho","dkappaa","dkappab","deta","dGcp","dGsp"]
            else:
                raise NotImplementedError(f"kernel type kltype={self._kltype} is not implemented for model type '{self._mdtype}'.")
        else:
            raise NotImplementedError(f"model type mdtype='{self._mdtype}' is not implemented. Supported types are 'iso' and 'dtti'.")
            
        return direc_list
    
    def grad_names(self,base=True):
        if self._mdtype == 'iso':
            grad_list = ['alpha_kernel','beta_kernel','rhop_kernel']
        else:
            grad_list = [f'c{i+1}{j+1}_kernel' for i in range(6) for j in range(i,6)]
            grad_list += ['rho_kernel']

        if not base:
            grad_list = self.direc_names()
            for i in range(len(grad_list)):
                grad_list[i] = grad_list[i][1:] + "_kernel"

        return grad_list
    
    def _cijkl2dtti(self,model,backward=False):
        # dispatch to auto-generated conversions (see auto_kergen.py)
        converter = MODEL_CONVERTERS.get(self._kltype)
        if converter is None:
            raise NotImplementedError(f"kernel type kltype={self._kltype} is not implemented for model type '{self._mdtype}'.")

        return converter(model,backward)

    def convert_model(self,model:np.ndarray,backward=False):
        """
        transform from base model to user defined model or vice versa

        Parameters
        -----------
        model: np.ndarray
            base model or user defined model
        backward: bool
            if True, convert from user defined model to base model

        Returns
        --------------
        model_new : np.ndarray
            user defined model or base model
        """
        if self._mdtype == "iso":
            model_new = model.copy()
            if self._kltype == 2:
                if not backward:
                    vp = model[0,...]
                    vs = model[1,...]
                    model_new[0,...] = vp / vs
                else:
                    vpvs = model[0,...]
                    vs = model[1,...]
                    model_new[0,...] = vpvs * vs 
        elif self._mdtype == "dtti":
            model_new = self._cijkl2dtti(model,backward)
        else:
            raise NotImplementedError(f"model type mdtype='{self._mdtype}' is not implemented. Supported types are 'iso' and 'dtti'.")
        
        return model_new
    
    def convert_to_visual(self,model:np.ndarray):
        """
        convert from base models to models for visualization

        Parameters
        -----------
        model: np.ndarray
            base model

        Returns
        --------------
        model_new : np.ndarray
            model for visualization
        plot_names: list[str]
            names for visualization
        """
        model_user = self.convert_model(model)
        plot_names = self.direc_names()
        for i in range(len(plot_names)):
            plot_names[i] = plot_names[i][1:]

        if self._mdtype == "dtti":
            if self._kltype == 1: # vp,vs,rho,gcp,gsp
                gcp = model_user[3,...]
                gsp = model_user[4,...]
                g0p = np.hypot(gsp,gcp)

                # phi are only in [-pi/2,pi/2]
                phi = 0.5 * np.arctan2(gsp,gcp)
                phi_deg = np.rad2deg(phi)
                phi_deg = (phi_deg + 90.) % 180. - 90.
                phi = np.deg2rad(phi_deg)
                
                # zero out phi when g0p < 1.0e-3
                idx = g0p < 1.0e-3
                phi[idx] = 0.

                # copy back to model_user
                model_user[3,...] = phi * 1.
                model_user[4,...] = g0p * 1.
                plot_names[3] = "phi"
                plot_names[4] = "G0"

            elif self._kltype == 2: # vph,vpv,vsh,vsv,rho,eta,gcp,gsp
                gcp = model_user[6,...]
                gsp = model_user[7,...]
                g0p = np.hypot(gsp,gcp)

                # phi are only in [-pi/2,pi/2]
                phi = 0.5 * np.arctan2(gsp,gcp)
                phi_deg = np.rad2deg(phi)
                phi_deg = (phi_deg + 90.) % 180. - 90.
                phi = np.deg2rad(phi_deg)

                # zero out phi when g0p < 1.0e-2
                idx = g0p < 1.0e-2
                phi[idx] = 0.

                # copy back to model_user
                model_user[6,...] = phi * 1.
                model_user[7,...] = g0p * 1.
                plot_names[6] = "phi"
                plot_names[7] = "G0"

            elif self._kltype == 3: # vp,vs,rho,(vph-vpv)/vpv,(vsh-vsv)/vsv, eta,gcp,gsp
                gcp = model_user[6,...]
                gsp = model_user[7,...]
                g0p = np.hypot(gsp,gcp)

                # phi are only in [-pi/2,pi/2]
                phi = 0.5 * np.arctan2(gsp,gcp)
                phi_deg = np.rad2deg(phi)
                phi_deg = (phi_deg + 90.) % 180. - 90.
                phi = np.deg2rad(phi_deg)

                # zero out phi when g0p < 1.0e-3
                idx = g0p < 1.0e-3
                phi[idx] = 0.
                
                # copy back to model_user
                model_user[6,...] = phi * 1.
                model_user[7,...] = g0p * 1.
                plot_names[6] = "phi"
                plot_names[7] = "G0"

            else:
                raise NotImplementedError(f"kernel type kltype={self._kltype} is not implemented for model type '{self._mdtype}'.")

        return model_user,plot_names

    def user2opt(self, md: np.ndarray):
        """
        Get the model used in gradient-based optimizers. Returns md for dimensionless parameters,
        and log(md) for others.

        Parameters
        -------------
        md: np.ndarray
            current user defined model
        
        Returns
        -------------
        md_used: np.ndarray
            model vector used in optimization
        """
        # initialize
        md_used = md * 1

        if self._mdtype == "iso":
            if self._kltype == 1: # vp,vs,rho
                md_used =  np.log(md)
            elif self._kltype == 2: #vp/vs,vs,rho
                md_used[1:,...] = np.log(md[1:,...])
        
        elif self._mdtype == "dtti":
            if self._kltype == 1: # vp,vs,rho,gcp,gsp
                md_used[:3,...] = np.log(md[:3,...])
            elif self._kltype == 2: # vph,vpv,vsh,vsv,rho,eta,gcp,gsp
                md_used[:5,...] = np.log(md[:5,...])
            elif self._kltype == 3: # vp,vs,rho,(vph-vpv)/vpv,(vsh-vsv)/vsv, eta,gcp,gsp
                md_used[:3,...] = np.log(md[:3,...])
            else:
                raise NotImplementedError(f"kernel type kltype={self._kltype} is not implemented for model type '{self._mdtype}'.")
        
        return md_used
    
    def model_update(self,md_usr:np.ndarray,direc:np.ndarray):
        """
        update model by given direction

        Parameters
        -------------
        md_usr: np.ndarray
            current user defined model
        direc: np.ndarray
            search direction in the optimization parameter space (see `user2opt`)
        
        Returns
        -------------
        md_update: np.ndarray
            updated user defined model
        """
        md_update = md_usr * 1.
        if self._mdtype == "iso":
            if self._kltype == 1: # vp,vs,rho
                md_update = md_usr * np.exp(direc)
            elif self._kltype == 2  : # vp/vs,vs,rho
                md_update[0,...] = md_usr[0,...] + direc[0,...]
                md_update[1:,...] = md_usr[1:,...] * np.exp(direc[1:,...])
        
        elif self._mdtype == "dtti":
            if self._kltype == 1: # vp,vs,rho,gcp,gsp
                md_update[3:,...] = md_usr[3:,...] + direc[3:,...]
                md_update[:3,...] = md_usr[:3,...] * np.exp(direc[:3,...])
            elif self._kltype == 2: # vph,vpv,vsh,vsv,rho,eta,gcp,gsp
                md_update[5:,...] = md_usr[5:,...] + direc[5:,...]
                md_update[:5,...] = md_usr[:5,...] * np.exp(direc[:5,...])
            elif self._kltype == 3: # vp,vs,rho,(vph-vpv)/vpv,(vsh-vsv)/vsv, eta,gcp,gsp
                md_update[3:,...] = md_usr[3:,...] + direc[3:,...]
                md_update[:3,...] = md_usr[:3,...] * np.exp(direc[:3,...])
            else:
                raise NotImplementedError(f"kernel type kltype={self._kltype} is not implemented for model type '{self._mdtype}'.")

        return md_update
    
    def _cijkl_kl2dtti(self,md_usr:np.ndarray,md_kl:np.ndarray):
        # dispatch to auto-generated conversions (see auto_kergen.py)
        converter = KERNEL_CONVERTERS.get(self._kltype)
        if converter is None:
            raise NotImplementedError(f"kernel type kltype={self._kltype} is not implemented for model type '{self._mdtype}'.")

        return converter(md_usr,md_kl)

    def convert_kl(self,md:np.ndarray,md_kl:np.ndarray):
        """
        convert base model to user defined model, and base kernels to
        kernels in the optimization parameter space

        Parameters
        -------------
        md: np.ndarray
            base model
        md_kl: np.ndarray
            kernels for base model

        Returns
        -------------
        md_usr: np.ndarray
            user defined model
        kl_opt: np.ndarray
            kernels in the optimization parameter space (see `user2opt`)
        """
        # first convert model to user defined model
        md_usr = self.convert_model(md)

        if self._mdtype == 'iso':
            kl_opt = md_kl.copy()
            if self._kltype == 2: # vpvs_vs_rho
                kl_opt[0,...] = md_kl[0,...] / md_usr[0,...]
                kl_opt[1,...] = md_kl[0,...] + md_kl[1,...] 
            pass

        elif self._mdtype == "dtti":
            kl_opt = self._cijkl_kl2dtti(md_usr,md_kl)
        else:
            raise NotImplementedError(f"model type mdtype='{self._mdtype}' is not implemented. Supported types are 'iso' and 'dtti'.")

        # mask part of the kernels
        kl_opt = self.mask_vector(kl_opt)


        return md_usr,kl_opt
    
    def mask_vector(self,kl:np.ndarray):
        """
        mask part of the kernels

        Parameters
        -------------
        kl: np.ndarray
            kernels in the optimization parameter space

        Returns
        -------------
        kl_masked: np.ndarray
            masked kernels in the optimization parameter space
        """
        kl_masked = kl * 1.
        kl_masked[self._mask_vars,...] = 0.

        return kl_masked