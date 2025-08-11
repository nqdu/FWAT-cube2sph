import os 
import shutil

from fwat.const import SOLVER,SRC_REC,PARAM_FILE
from fwat.system.specfem import get_param,change_parfile
from fwat.system.system import sys_remove

def get_noise_mc_event_list(evtid:str,verbose = True):
    # check  required CC_COMPS
    import yaml 

    with open(PARAM_FILE,"r") as fio:
        cc_comps = yaml.safe_load(fio)['measure']['noise']['CC_COMPS']

    if verbose: print("")
    evtid_set = set()
    scomps = sorted([cc[0] for cc in cc_comps])
    for ch in ['E','N','Z']:
        name = f"{SRC_REC}/STATIONS_{evtid}_{ch}"
        if ch in scomps and os.path.exists(name):
            evtid_set.add(f"{evtid}_{ch}")

        # check if station_files is missing
        if ch in scomps and (not os.path.exists(name)) and verbose:
            print(f"multi-channel noise evtid = {evtid}:")
            print(f"source channel {ch} is required but station file {name} is missing")
            print(f"skip source channel {ch}")

    for ch in ['R','T']:
        name1 = f"{SRC_REC}/STATIONS_{evtid}_{ch}"
        if ch in scomps and os.path.exists(name1):
            evtid_set.add(f"{evtid}_E")
            evtid_set.add(f"{evtid}_N")

        if ch in scomps and (not os.path.exists(name1)) and verbose:
            print(f"multi-channel noise evtid = {evtid}:")
            print(f"source channel {ch} is required but station file {name1} is missing")
            print(f"skip source channel {ch}")

    # get unique
    evtid_list = sorted(list(evtid_set))
    
    return evtid_list


class FwatSubmitor:
    def __init__(self,meatype:str,iter0:int,evtid:str,run_opt:int):        
        # load submit parameters
        import yaml
        with open(f"{PARAM_FILE}","r") as f:
            param = yaml.safe_load(f)
            self.pdict = param['simulation']
            self.mdict = param['measure'][meatype]

        # get GPU_MODE in DATA/Par_file.{meatype}
        from fwat.system.specfem import get_param
        mode = get_param(f'DATA/Par_file.{meatype}','GPU_MODE').lower()
        if mode == '.false.':
            self.pdict['GPU_MODE'] = False 
        else:
            self.pdict['GPU_MODE'] = True

        # get path
        self.SOLVER = SOLVER
        self.SRC_REC = SRC_REC

        # current working directory
        self.cwd = os.getcwd()

        # check measure type
        assert(meatype in ['tele','sks','rf','noise'])

        # path to simulation directory
        mdir = "M%02d" %(iter0)
        if run_opt == 2:
            mdir = mdir + ".ls"
        self.evtid = evtid
        self.meatype = meatype
        self.mod = mdir
        self.run_opt = run_opt

        # nprocs and path
        self.nprocs:int = get_param(f'./DATA/Par_file.{self.meatype}','NPROC')

    def _prepare_model_and_parfile(self,evtid:str):
        # create directories
        syndir:str = f"{self.SOLVER}/{self.mod}/{evtid}/"
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
        if self.run_opt == 1:
            direc = f"{self.cwd}/{LOCAL_PATH}"
        else:
            direc = f"{self.cwd}/optimize/MODEL_{self.mod}/"

        for f in os.listdir(direc):
            os.symlink(f"{direc}/{f}",f"{syndir}/{LOCAL_PATH}/{f}")

    def _filter_station_files(self,file_list):
        with open(file_list[0]) as f:
            all_lines = f.readlines()
        for fname in file_list[1:]:
            with open(fname) as f:
                all_lines.extend(f.readlines())

        unique_sorted = sorted(set(line.strip() for line in all_lines))

        return unique_sorted

    def _get_simulation_list(self,verbose=True):
        # get how many simulations required
        evtid_list = []
        stations_list = []
        if self.meatype != 'noise':
            evtid_list.append(self.evtid)
            name = f"{self.SRC_REC}/STATIONS_{self.evtid}"
            strout = self._filter_station_files([name])
            stations_list.append(strout)
        else:
            evtid_list = get_noise_mc_event_list(self.evtid,verbose=verbose)
            if len(evtid_list) == 0:
                print(f"for noise simulation evtid = {self.evtid}, no valid source components find!")
                exit(1)
            
            # get stations_list
            import glob 
            name = f"{self.SRC_REC}/STATIONS_{self.evtid}_[ENZRT]"
            names = glob.glob(name)
            strout = self._filter_station_files(names)
            for _ in range(len(evtid_list)):
                stations_list.append(strout)
        
        return evtid_list,stations_list

    def prepare_fwd(self):
        import glob 

        evtid_list,stations_list = self._get_simulation_list(True)

        # loop each evtid list
        nevts = len(evtid_list)
        for ie in range(nevts):
            evtid = evtid_list[ie]
            print(f"prepare parameters for {self.meatype}, evtid = {evtid}")
            self._prepare_model_and_parfile(evtid_list[ie])

            # prepare stations
            syndir = f"{self.SOLVER}/{self.mod}/{evtid}"
            with open(f"{syndir}/DATA/STATIONS","w") as f:
                for line in stations_list[ie]:
                    f.write(line + "\n")
            shutil.copy2(f"{syndir}/DATA/STATIONS",
                         f"{syndir}/DATA/STATIONS_ADJOINT")

            # prepare source
            if self.meatype in ['tele','sks','rf']:
                # filename
                filename = f"{syndir}/DATA/Par_file"

                # change parfile
                change_parfile(filename,
                                COUPLE_WITH_INJECTION_TECHNIQUE=".true.",
                                INJECTION_TECHNIQUE_TYPE=4)

                # create an empty forcesolution file
                with open(f"{syndir}/DATA/FORCESOLUTION","w") as fio:
                    pass

                # link axisem field
                for f in glob.glob("f{syndir}/DATABASES_MPI/*wavefield_discontinuity.bin"):
                    sys_remove(f)
                for f in os.listdir(f"{self.cwd}/DATA/axisem/{evtid}/"):
                    os.symlink(f"{self.cwd}/DATA/axisem/{evtid}/{f}",
                               f"{syndir}/DATABASES_MPI/{f}" )
            
            elif self.meatype == "noise":
                filename = f"{syndir}/DATA/Par_file"
                change_parfile(filename,
                                COUPLE_WITH_INJECTION_TECHNIQUE=".false.",
                                INJECTION_TECHNIQUE_TYPE=4)

                # copy FORCESOLUTION FILE
                filename = f"{self.SRC_REC}/FORCESOLUTION_{evtid}"
                shutil.copy2(filename,f"{syndir}/DATA/FORCESOLUTION")

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

            # if run_opt == 1 diasble SUBSAMPLE
            if self.run_opt == 1:
                 SUBSAMPLE = ".false."

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
            
        return evtid_list

    def prepare_adj(self):
        evtid_list,_ = self._get_simulation_list(False)

        # loop each evtid list
        nevts = len(evtid_list)
        for ie in range(nevts):
            evtid = evtid_list[ie]
            syndir = f"{self.SOLVER}/{self.mod}/{evtid}"

            # change pars
            GPU_MODE = '.true.'
            if not self.pdict['GPU_MODE']:
                GPU_MODE = '.false.'
            
            # subsample
            SUBSAMPLE = ".false."
            COUPLE_WITH_INJECTION = ".false."
            SAVE_FORWARD=".false."
            if self.pdict['DUMP_WAVEFIELDS']:
                SUBSAMPLE = ".true."
            else:
                if self.meatype in ['sks','tele','rf']:
                    COUPLE_WITH_INJECTION = ".true."
            
            filename = f"{syndir}/DATA/Par_file"
            change_parfile(filename,
                        COUPLE_WITH_INJECTION_TECHNIQUE=COUPLE_WITH_INJECTION,
                        SUBSAMPLE_FORWARD_WAVEFIELD=SUBSAMPLE,
                        SAVE_FORWARD=SAVE_FORWARD,
                        GPU_MODE=GPU_MODE,
                        SIMULATION_TYPE=3,
                        APPROXIMATE_HESS_KL=".true.")

        return evtid_list

    def _get_timestamp():
        from datetime import datetime

        # Get the current date and time
        current_datetime = datetime.now()

        # You can also format the output using strftime()
        formatted_datetime = current_datetime.strftime("%Y-%m-%d %H:%M:%S")
        
        return formatted_datetime


def prepare_fwd(args):
    if len(args) != 4:
        print("Usage: prepare_fwd fwat meatype,iter0,evtid,run_opt")
    
    meatype:str = args[0]
    iter0:int = int(args[1])
    evtid:str = args[2]
    run_opt:int = int(args[3])

    op = FwatSubmitor(meatype,iter0,evtid,run_opt)
    return op.prepare_fwd()

def prepare_adj(args):
    if len(args) != 4:
        print("Usage: prepare_adj fwat meatype,iter0,evtid,run_opt")

    
    meatype:str = args[0]
    iter0:int = int(args[1])
    evtid:str = args[2]
    run_opt:int = int(args[3])

    op = FwatSubmitor(meatype,iter0,evtid,run_opt)
    return op.prepare_adj()