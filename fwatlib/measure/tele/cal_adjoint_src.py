import numpy as np
import sys
from obspy.io.sac import SACTrace
from scipy.signal import correlate

def main():
    if len(sys.argv) != 4:
        print("Usage: python cal_adjoint_src.py adjfile stf_sac avgamp")
        exit(1)
    
    # get parma
    adjfile = sys.argv[1]
    stf_sac = sys.argv[2]
    avgamp = float(sys.argv[3])

    # read stf 
    tr = SACTrace.read(stf_sac)
    dt = tr.delta

    # read data
    d = np.loadtxt(adjfile)

    # compute adjoint source
    x  = d[:,1] * 1.
    x /= avgamp
    tmp = correlate(x,tr.data,'same') * dt 
    if np.max(np.abs(tmp)) > 0.:
        x = tmp / np.max(np.abs(tr.data))

    # save data
    d[:,1] = x
    np.savetxt(adjfile,d)


if __name__ == "__main__":
    main()