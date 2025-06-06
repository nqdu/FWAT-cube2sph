import os 

class FwatSubmitor:
    def __init__(self,module_env_file:str):
        # load module env
        with open(module_env_file,"r") as fio:
            self._env_cmd = fio.readlines()
        
        # load submit parameters
        import yaml
        with open("fwat_params/FWAT.PAR.yaml","r") as f:
            self.pdict = yaml.safe_load(f)['submit']
    
    def prepare_files(self,meatype:str,iter0:int,evtid:int,run_opt:int):
        # path to simulation directory
        mdir = "M%02d" %(iter0)
        if self.run_opt == 2:
            mdir = mdir + ".ls"
        self.mod = mdir
        self.syndir:str = f"{self.SOLVER}/{mdir}/{evtid}/"

        # create directories
        os.makedirs(f"{self.syndir}")
        os.makedirs(f"{self.syndir}/DATA")
        os.makedirs(f"{self.syndir}/OUTPUT_FILES")
        os.makedirs(f"{self.syndir}")

        