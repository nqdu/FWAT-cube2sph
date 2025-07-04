import numpy as np 

def _Index(m,n):
    assert(m<=n)
    idx = m * 6 + n - (m * (m + 1)) // 2

    return idx

class FwatModel:
    def initialize(self,mdtype='iso',kltype=2) -> None:
        """
        initialize FwatModel

        Parameters
        ----------

        mdtype: str
            currently support [iso,dtti]
        kltype: int
            iso configuration
                0: kappa,mu,rho
                1: vp,vs,rho
                2: vp/vs,vs,rho
            dtti configuration:
                0: c11-c66,rho
                1: vp,vs,rho,gcp,gsp
                2. vph,vpv,vsh,vsv,rho,eta,gcp,gsp
        """
        self._mdtype = mdtype
        self._kltype =kltype

        if mdtype not in ['iso','dtti']:
            print(f"not implemented for modeltype = {mdtype}")
            exit(1)
    
    def __init__(self,filename='fwat_params/FWAT.PAR.yaml',mdtype:str = None,kltype:int = None) -> None:
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
        else:
            assert(mdtype is not None)
            assert(kltype is not None)
    
        self.initialize(mdtype,kltype)

    def get_model_names(self):
        if self._mdtype == 'iso':
            mod_list = ['vp','vs','rho']
        else:
            mod_list = [f'c{i+1}{j+1}' for i in range(6) for j in range(i,6)]
            mod_list += ['rho']
        
        return mod_list
    
    def get_direc_names(self):
        if self._mdtype == 'iso':
            direc_list = ['dalpha','dbeta','drho']
            #grad_list = ['alpha_kernel_smooth','beta_kernel_smooth','rhop_kernel_smooth']
            if self._kltype == 2:
                direc_list[0] = 'dvpvs'
        elif self._mdtype == "dtti":
            direc_list = [f'dc{i+1}{j+1}' for i in range(6) for j in range(i,6)]
            direc_list += ['drho']
            if self._kltype == 1:
                # vp,vs,rho,gc_nodim,gs_nodim  Zhu et al 2015, GJI, (25,26) 
                direc_list = ["dalpha","dbeta","drho","dGcp","dGsp"]
            elif self._kltype == 2:
                direc_list = ["dalphah","dbetav","dalphah","dbetav","drho","deta","dGcp","dGsp"]

            
        return direc_list
    
    def get_grad_names(self,base=True):
        if self._mdtype == 'iso':
            grad_list = ['alpha_kernel','beta_kernel','rhop_kernel']
        else:
            grad_list = [f'c{i+1}{j+1}_kernel' for i in range(6) for j in range(i,6)]
            grad_list += ['rho_kernel']
            if not base:
                grad_list = self.get_direc_names()
                for i in range(len(grad_list)):
                    grad_list[i] = grad_list[i][1:] + "_kernel"

        return grad_list
    
    def _cijkl2dtti(self,model,backward=False):
        # model now is c11, c12 -- c66, rho with shape(22）
        if self._kltype == 1: # vp,vs,rho,gc gs
            if not backward:
                rho = model[-1,:] * 1.

                # model with shape(22,...)
                new_shape = (5,) + model.shape[1:]
                model_new = np.zeros(new_shape)

                # get model
                vp = np.sqrt(model[_Index(0,0),...] / rho)
                vs = np.sqrt(model[_Index(5,5),...] / rho)
                gc = 0.5 * (model[_Index(4,4),...] - model[_Index(3,3),...])
                gs = - model[_Index(3,4),...]

                # vp,vs,rho,gc_nodim,gs_nodim
                model_new[0,...] = vp * 1.
                model_new[1,...] = vs * 1.
                model_new[2,...] = rho 
                model_new[3,...] = gc  / (rho * vs**2)
                model_new[4,...] = gs  / (rho * vs**2)
            else:
                # model with shape(5,...)
                new_shape = (22,) + model.shape[1:]
                model_new = np.zeros(new_shape)

                vp  = model[0,...] * 1.
                vs  = model[1,...] * 1.
                rho = model[2,...] * 1.
                gc = (rho * vs**2) * model[3,...]
                gs = (rho * vs**2) * model[4,...]

                # convert to tti
                A = vp**2 * rho 
                L = vs**2 * rho 
                N = vs**2 * rho 
                C = vp**2 * rho 
                F = A - 2 * L 

                model_new[-1,...] = rho 
                model_new[_Index(0,0),...] = A
                model_new[_Index(0,1),...] = A - 2 * N  
                model_new[_Index(0,2),...] = F
                model_new[_Index(1,1),...] = A
                model_new[_Index(1,2),...] = F
                model_new[_Index(2,2),...] = C
                model_new[_Index(3,3),...] = L - gc
                model_new[_Index(3,4),...] = -gs
                model_new[_Index(4,4),...] = L + gc
                model_new[_Index(5,5),...] = N

        elif self._kltype == 2: # vph,vpv,vsh,vsv,rho,eta,gcp,gsp
            if not backward:
                rho = model[-1,:] * 1.
                new_shape = (8,) + model.shape[1:]
                model_new = np.zeros(new_shape)

                # get parameters
                A = model[_Index(0,0),...]
                C = model[_Index(1,1),...]
                L = 0.5 * (model[_Index(3,3),...] + model[_Index(4,4),...])
                N = model[_Index(5,5),...]
                F = model[_Index(0,2),...]
                eta = F / (A - 2 * L)
                gcp = (model[_Index(4,4),...] - L) / L
                gsp = -model[_Index(3,4),...] / L 

                # save to model_new
                model_new[0,...] = np.sqrt(A / rho) # vph
                model_new[1,...] = np.sqrt(C / rho) # vpv
                model_new[2,...] = np.sqrt(N / rho) # vsh
                model_new[3,...] = np.sqrt(L / rho) # vsv
                model_new[4,...] = rho * 1.
                model_new[5,...] = eta * 1. 
                model_new[6,...] = gcp 
                model_new[7,...] = gsp 

                pass
            else:
                new_shape = (22,) + model.shape[1:]
                model_new = np.zeros(new_shape)

                # get parameters
                rho = model[4,...]
                vph = model[0,...]; A = vph**2 * rho 
                vpv = model[1,...]; C = vpv**2 * rho 
                vsh = model[2,...]; N = vsh**2 * rho 
                vsv = model[3,...]; L = vsv**2 * rho
                eta = model[5,...]
                gc = model[6,...] * L
                gs = model[7,...]* L 
                F = eta * (A - 2 * L)

                # copy to cijkl
                model_new[_Index(0,0),...] = A
                model_new[_Index(1,1),...] = A
                model_new[_Index(2,2),...] = C 
                model_new[_Index(0,1),...] = A - 2 * N 
                model_new[_Index(0,2),...] = F
                model_new[_Index(1,2),...] = F
                model_new[_Index(3,3),...] = L - gc 
                model_new[_Index(4,4),...] = L + gc 
                model_new[_Index(5,5),...] = N 
                model_new[_Index(3,4),...] = -gs
                model_new[-1,...] = rho

        return model_new
    
    def convert_md(self,model:np.ndarray,backward=False):
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
        
        return model_new
    
    def convert_md_visual(self,model:np.ndarray):
        """
        convert from base models to modesl for visualization

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
        model_new = self.convert_md(model)
        plot_names = self.get_direc_names()
        for i in range(len(plot_names)):
            plot_names[i] = plot_names[i][1:]

        if self._mdtype == "dtti":
            if self._kltype == 1: # vp,vs,rho,gcp,gsp
                gcp = model_new[3,...]
                gsp = model_new[4,...]
                phi = 0.5 * np.arctan2(gsp,gcp)
                g0p = np.hypot(gsp,gcp)

                # copy back to model_new
                model_new[3,...] = np.float32(phi)
                model_new[4,...] = np.float32(g0p)
                plot_names[3] = "phi"
                plot_names[4] = "G0"

            elif self._kltype == 2: # vph,vpv,vsh,vsv,rho,eta,gcp,gsp
                gcp = model_new[6,...]
                gsp = model_new[7,...]
                phi = 0.5 * np.arctan2(gsp,gcp)
                g0p = np.hypot(gsp,gcp)

                # copy back to model_new
                model_new[5,...] = np.float32(phi)
                model_new[6,...] = np.float32(g0p)
                plot_names[5] = "phi"
                plot_names[6] = "G0"

        return model_new,plot_names
    
    def get_used_model(self,md:np.ndarray):
        """
        get model used in gradient based optimizer, return md for dimensionless parameter,
        and log(md) for others
        """
        # initilize 
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
        
        return md_used
    
    def model_update(self,md:np.ndarray,direc:np.ndarray):
        """
        update model by given direction
        """
        md_update = md * 1.
        if self._mdtype == "iso":
            if self._kltype == 1: # vp,vs,rho
                md_update = md * np.exp(direc)
            elif self._kltype == 2: # vp/vs,vs,rho
                md_update[0,...] = md[0,...] + direc[0,...]
                md_update[1:,...] = md[1:,...] * np.exp(direc[1:,...])
        
        elif self._mdtype == "dtti":
            if self._kltype == 1: # vp,vs,rho,gcp,gsp
                md_update[3:,...] = md[3:,...] + direc[3:,...]
                md_update[:3,...] = md[:3,...] * np.exp(direc[:3,...])
            elif self._kltype == 2: # vph,vpv,vsh,vsv,rho,eta,gcp,gsp
                md_update[5:,...] = md[5:,...] + direc[5:,...]
                md_update[:5,...] = md[:5,...] * np.exp(direc[:5,...])

        return md_update
    
    def  _cijkl_kl2dtti(self,md_new:np.ndarray,md_kl:np.ndarray):
        if self._kltype == 1:
            vp  = md_new[0,...] * 1.
            vs  = md_new[1,...] * 1.
            rho = md_new[2,...] * 1.
            gcp = md_new[3,...]
            gsp = md_new[4,...]
            md_newkl = md_new * 0

            # auto generated by sympy
            md_newkl[0,...] = md_kl[0,...]*(2*rho*vp)+md_kl[1,...]*(2*rho*vp)+md_kl[2,...] \
                *(2*rho*vp)+md_kl[6,...]*(2*rho*vp)+md_kl[7,...]*(2*rho*vp)+ \
                md_kl[11,...]*(2*rho*vp)
            md_newkl[1,...] = md_kl[1,...]*(-4*rho*vs)+md_kl[2,...]*(-4*rho*vs)+ \
                md_kl[7,...]*(-4*rho*vs)+md_kl[15,...]*(-2*gcp*rho*vs + 2* \
                rho*vs)+md_kl[16,...]*(-2*gsp*rho*vs)+md_kl[18,...]*(2*gcp* \
                rho*vs + 2*rho*vs)+md_kl[20,...]*(2*rho*vs)
            md_newkl[2,...] = md_kl[-1,...]* (1)+md_kl[0,...]*(vp**2)+md_kl[1,...]*(vp**2 \
                - 2*vs**2)+md_kl[2,...]*(vp**2 - 2*vs**2)+md_kl[6,...]*(vp** \
                2)+md_kl[7,...]*(vp**2 - 2*vs**2)+md_kl[11,...]*(vp**2)+ \
                md_kl[15,...]*(-gcp*vs**2 + vs**2)+md_kl[16,...]*(-gsp*vs** \
                2)+md_kl[18,...]*(gcp*vs**2 + vs**2)+md_kl[20,...]*(vs**2)
            md_newkl[3,...] = md_kl[15,...]*(-rho*vs**2)+md_kl[18,...]*(rho*vs**2)
            md_newkl[4,...] = md_kl[16,...]*(-rho*vs**2)

            # get relative change if required
            md_newkl[0:3,...] *= md_new[0:3,...]
            
        elif self._kltype == 2:
            vph = md_new[0,...] * 1.
            vpv = md_new[1,...] * 1.
            vsh  = md_new[2,...] * 1.
            vsv  = md_new[3,...] * 1.
            rho = md_new[4,...] * 1.
            eta = md_new[5,...] * 1.
            gcp = md_new[6,...]
            gsp = md_new[7,...]

            # auto generated by sympy
            md_newkl[0,...] = md_kl[0,...]*(2*rho*vph)+md_kl[1,...]*(2*rho*vph)+ \
                md_kl[2,...]*(2*eta*rho*vph)+md_kl[6,...]*(2*rho*vph)+ \
                md_kl[7,...]*(2*eta*rho*vph)
            md_newkl[1,...] = md_kl[11,...]*(2*rho*vpv)
            md_newkl[2,...] = md_kl[1,...]*(-4*rho*vsh)+md_kl[20,...]*(2*rho*vsh)
            md_newkl[3,...] = md_kl[2,...]*(-4*eta*rho*vsv)+md_kl[7,...]*(-4*eta*rho*vsv)+ \
                md_kl[15,...]*(-2*gcp*rho*vsv + 2*rho*vsv)+md_kl[16,...]*(-2 \
                *gsp*rho*vsv)+md_kl[18,...]*(2*gcp*rho*vsv + 2*rho*vsv)
            md_newkl[4,...] = md_kl[-1,...]* (1)+md_kl[0,...]*(vph**2)+md_kl[1,...]*(vph** \
                2 - 2*vsh**2)+md_kl[2,...]*(eta*(vph**2 - 2*vsv**2))+ \
                md_kl[6,...]*(vph**2)+md_kl[7,...]*(eta*(vph**2 - 2*vsv**2)) \
                +md_kl[11,...]*(vpv**2)+md_kl[15,...]*(-gcp*vsv**2 + vsv**2) \
                +md_kl[16,...]*(-gsp*vsv**2)+md_kl[18,...]*(gcp*vsv**2 + vsv \
                **2)+md_kl[20,...]*(vsh**2)
            md_newkl[5,...] = md_kl[2,...]*(rho*vph**2 - 2*rho*vsv**2)+md_kl[7,...]*(rho* \
                vph**2 - 2*rho*vsv**2)
            md_newkl[6,...] = md_kl[15,...]*(-rho*vsv**2)+md_kl[18,...]*(rho*vsv**2)
            md_newkl[7,...] = md_kl[16,...]*(-rho*vsv**2)

            # get relative change if required
            md_newkl[0:5,...] *= md_new[0:5,...]

        return md_newkl
    
    def convert_kl(self,md:np.ndarray,md_kl:np.ndarray):
        """
        convert original kernels/models to parameter set used

        Parameters
        -------------
        md: np.ndarray
            base model
        md_kl: np.ndarray
            kernels for base model

        Returns
        -------------
        md_new: np.ndarray
            user defined model
        md_newkl: np.ndarray
            kernels for user defined model
        """
        # first convert model to base one
        md_new = self.convert_md(md)

        if self._mdtype == 'iso':
            md_newkl = md_kl.copy()
            if self._kltype == 2: # vpvs_vs_rho
                md_newkl[0,...] = md_kl[0,...] / md_new[0,...]
                md_newkl[1,...] = md_kl[0,...] + md_kl[1,...] 
            pass

        elif self._mdtype == "dtti":
            md_newkl = self._cijkl_kl2dtti(md_new,md_kl)


        # mask part of the kernels
        md_newkl[self._mask_vars,...] = 0.

        return md_new,md_newkl