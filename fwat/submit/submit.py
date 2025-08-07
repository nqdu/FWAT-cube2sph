import os 
import shutil

from fwat.const import SOLVER,SRC_REC,PARAM_FILE
from fwat.system.specfem import get_param,change_parfile
from fwat.system.system import sys_remove

class FwatSubmitor:
    def __init__(self,meatype:str,iter0:int,evtid:str,run_opt:int):        
        # load submit parameters
        import yaml
        with open(f"{PARAM_FILE}","r") as f:
            self.pdict = yaml.safe_load(f)['simulation']

        # get path
        self.SOLVER = SOLVER
        self.SRC_REC = SRC_REC

        # current working directory
        self.cwd = os.getcwd()

        # path to simulation directory
        mdir = "M%02d" %(iter0)
        if run_opt == 2:
            mdir = mdir + ".ls"
        syndir:str = f"{self.SOLVER}/{mdir}/{evtid}/"
        self.syndir = syndir
        self.evtid = evtid
        self.meatype = meatype
        self.mod = mdir

        # check if evtid is 

        # nprocs
        self.nprocs = get_param(f'./DATA/Par_file.{self.meatype}','NPROC')

    def _sanity_check(self):
        if self.meatype == "noise":
            pass
        pass 
    
    def prepare_files_fwd(self):        
        # create directories
        syndir = self.syndir
        LOCAL_PATH = './DATABASES_MPI/'
        os.makedirs(f"{syndir}",exist_ok=True)
        os.makedirs(f"{syndir}/DATA",exist_ok=True)
        os.makedirs(f"{syndir}/OUTPUT_FILES",exist_ok=True)
        os.makedirs(f"{syndir}/{LOCAL_PATH}",exist_ok=True)

        # copy file to syndir
        shutil.copy2(f'DATA/Par_file.{self.meatype}',f'{syndir}/DATA/Par_file')
        for f in os.listdir("OUTPUT_FILES/"):
            if '.h' in f:
                shutil.copy2(f"OUTPUT_FILES/{f}",f'{syndir}/OUTPUT_FILES/{f}')

        # create softlink for mesh files
        if os.path.exists(f"{syndir}/DATA/meshfem3D_files"):
            sys_remove(f"{syndir}/DATA/meshfem3D_files")
        os.makedirs(f"{syndir}/DATA/meshfem3D_files")
        for f in os.listdir(f"DATA/meshfem3D_files/"):
            name = f"{syndir}/DATA/meshfem3D_files/{f}"
            os.symlink(f"{self.cwd}/DATA/meshfem3D_files/{f}", name)

        # change local path
        filename = f"{syndir}/DATA/Par_file"
        change_parfile(filename,LOCAL_PATH=LOCAL_PATH)
        filename = f"{syndir}/DATA/meshfem3D_files/Mesh_Par_file"
        change_parfile(filename,LOCAL_PATH=LOCAL_PATH)

        # create model link
        if os.path.exists(f"{syndir}/{LOCAL_PATH}"):
            sys_remove(f"{syndir}/{LOCAL_PATH}")
        os.makedirs(f"{syndir}/{LOCAL_PATH}")
        for f in os.listdir(f"{self.cwd}/optimize/MODEL_{self.mod}/"):
            os.symlink(f"{self.cwd}/optimize/MODEL_{self.mod}/{f}",f"{syndir}/{LOCAL_PATH}/{f}")

        # copy stations
        evtid = self.evtid
        shutil.copy2(f"{self.SRC_REC}/STATIONS_{evtid}",
                     f"{syndir}/DATA/STATIONS")
        shutil.copy2(f"{self.SRC_REC}/STATIONS_{evtid}",
                     f"{syndir}/DATA/STATIONS_ADJOINT")
        
        # change forward parameters
        meatype = self.meatype
        if meatype in ['tele','sks','rf']:
            self._prepare_tele()
        elif meatype == "noise":
            self._prepare_noise()
            
        else:
            print("not implemented!")
            exit(1)

        # get simulation parameters
        GPU_MODE = '.true.'
        if not self.pdict['GPU_MODE']:
            GPU_MODE = '.false.'
        SUBSAMPLE = ".false."
        SAVE_FORWARD=".true."
        if self.pdict['DUMP_WAVEFIELDS']:
            SUBSAMPLE = ".true."
            SAVE_FORWARD=".false."

        # forward simulation
        filename = f"{syndir}/DATA/Par_file"
        NSTEP = get_param(filename,"NSTEP")
        change_parfile(filename,
                       SUBSAMPLE_FORWARD_WAVEFIELD=SUBSAMPLE,
                       SAVE_FORWARD=SAVE_FORWARD,
                       GPU_MODE=GPU_MODE,
                       SIMULATION_TYPE=1,
                       APPROXIMATE_HESS_KL=".false.",
                       WRITE_SEISMOGRAMS_BY_MASTER=".TRUE.",
                       SAVE_ALL_SEISMOS_IN_ONE_FILE=".true.",
                       NTSTEP_BETWEEN_OUTPUT_SEISMOS=NSTEP)

    def _get_timestamp():
        from datetime import datetime

        # Get the current date and time
        current_datetime = datetime.now()

        # You can also format the output using strftime()
        formatted_datetime = current_datetime.strftime("%Y-%m-%d %H:%M:%S")
        
        return formatted_datetime

    def prepare_files_adj(self):
        # get simulation parameters
        GPU_MODE = '.true.'
        if not self.pdict['GPU_MODE']:
            GPU_MODE = '.false.'
        
        # subsample
        SUBSAMPLE = ".false."
        COUPLE_WITH_INJECTION = ".false."
        SAVE_FORWARD=".false."
        if self.pdict['DUMP_WAVEFIELDS']:
            SUBSAMPLE = ".true."
        
        filename = f"{self.syndir}/DATA/Par_file"
        change_parfile(filename,
                       COUPLE_WITH_INJECTION_TECHNIQUE=COUPLE_WITH_INJECTION,
                       SUBSAMPLE_FORWARD_WAVEFIELD=SUBSAMPLE,
                       SAVE_FORWARD=SAVE_FORWARD,
                       GPU_MODE=GPU_MODE,
                       SIMULATION_TYPE=3,
                       APPROXIMATE_HESS_KL=".true.")

    def _prepare_tele(self):
        # filename
        syndir = self.syndir
        filename = f"{self.syndir}/DATA/Par_file"

        # change parfile
        change_parfile(filename,
                        COUPLE_WITH_INJECTION_TECHNIQUE=".true.",
                        INJECTION_TECHNIQUE_TYPE=4    )

        # create an empty forcesolution file
        with open(f"{syndir}/DATA/FORCESOLUTION","w") as fio:
            pass

        # link axisem field
        for f in os.listdir(f"{syndir}/DATABASES_MPI/"):
            if 'wavefield_discontinuity' in f :
                sys_remove(f"{syndir}/DATABASES_MPI/{f}")
                os.symlink(f"{self.cwd}/DATA/axisem/{self.evtid}/{f}", 
                          f"{syndir}/DATABASES_MPI/{f}")
                
    def _prepare_noise(self):
        # change parfile
        filename = f"{self.syndir}/DATA/Par_file"
        change_parfile(filename,
                        COUPLE_WITH_INJECTION_TECHNIQUE=".false.",
                        INJECTION_TECHNIQUE_TYPE=4)

        # copy FORCESOLUTION FILE
        filename = f"{self.SRC_REC}/FORCESOLUTION_{self.evtid}"
        shutil.copy2(filename,f"{self.syndir}/DATA/FORCESOLUTION")


def prepare_fwd(args):
    if len(args) != 4:
        print("Usage: prepare_fwd fwat meatype,iter0,evtid,run_opt")
    
    meatype:str = args[0]
    iter0:int = int(args[1])
    evtid:str = args[2]
    run_opt:int = int(args[3])

    op = FwatSubmitor(meatype,iter0,evtid,run_opt)
    op.prepare_files_fwd()

def prepare_adj(args):
    if len(args) != 4:
        print("Usage: prepare_adj fwat meatype,iter0,evtid,run_opt")
    
    meatype:str = args[0]
    iter0:int = int(args[1])
    evtid:str = args[2]
    run_opt:int = int(args[3])

    op = FwatSubmitor(meatype,iter0,evtid,run_opt)
    op.prepare_files_adj()