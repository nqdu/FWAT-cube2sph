import numpy as np
from scipy.interpolate import RegularGridInterpolator
import sys 
from scipy.io import FortranFile
from tools import *

def coswin(n, opt):
    nn = 2 * n + 1
    win = np.zeros((nn))

    for i in range(nn):
        win[i] = 0.5 - 0.5 * np.cos(2 * np.pi * i / (nn-1))
    
    if opt == 1:
        winx = win[:n]
    else:
        winx = win[n+1:nn]
    
    return winx

def create_taper1d(taph1,taph2,nx,hx):
    nh1 = int(np.floor(taph1/hx))
    nh2 = int(np.floor(taph2/hx))
    x1 = nh1
    x2 = nh1 + nh2
    x3 = nx - nh1 - nh2
    x4 = nx - nh1
    xtaper = np.zeros((nx))
    xtaper[x1:x1+nh2] = coswin(nh2,1)
    xtaper[x3:x4] = coswin(nh2,0)
    xtaper[x2:x3] = 1.

    return xtaper 

def main():
    if len(sys.argv) !=3 :
        print("need 6 parameters: iter nprocs")
        print("example: python taper_kernel.py 0 160")
        exit(1)

    iter = int(sys.argv[1])
    nprocs = int(sys.argv[2])

    # taper params
    f = open("fwat_params/FWAT.PAR","r")
    lines = f.readlines()
    f.close()
    taph1 = 0.; taph2 = 0.; tapz1 = 0.; tapz2 = 0.
    for line in lines:
        if 'TAPER_H_SUPPRESS' in line:
            taph1 = float(line.split(':')[-1])
        if 'TAPER_H_BUFFER' in line:
            taph2 = float(line.split(':')[-1])
        if 'TAPER_V_SUPPRESS' in line:
            tapz1 = float(line.split(':')[-1])
        if 'TAPER_V_BUFFER' in line:
            tapz2 = float(line.split(':')[-1])

    # read model info
    f = open("OUTPUT_FILES/output_solver.txt","r")
    lines = f.readlines()
    f.close()
    elem_min = 0.; xmin = 0.; xmax = 0
    ymin = 0.; ymax = 0.; zmin = 0; zmax = 0
    for line in lines:
        if 'Min element size' in line:
            elem_min = np.float32(line.split('=')[-1])
        if 'Xmin' in line:
            info =  line.split('=')[-1].split()
            xmin = float(info[0]); xmax = float(info[1])
        if 'Ymin' in line:
            info =  line.split('=')[-1].split()
            ymin = float(info[0]); ymax = float(info[1])
        if 'Zmin' in line:
            info =  line.split('=')[-1].split()
            zmin = float(info[0]); zmax = float(info[1])

    print(xmin,xmax,ymin,ymax,zmin,zmax,elem_min)

    # create arrays
    hx = elem_min * 0.5
    hy = hx; hz = hx
    nx = int(np.floor((xmax-xmin)/hx+1))
    ny = int(np.floor((ymax-ymin)/hy+1))
    nz = int(np.floor((zmax-zmin)/hz+1))
    x = xmin + np.arange(nx) * hx 
    y = ymin + np.arange(ny) * hy
    z = zmin + np.arange(nz) * hz
    print(nx,ny,nz)

    # create taper
    xtaper = create_taper1d(taph1,taph2,nx,hx)
    ytaper = create_taper1d(taph1,taph2,ny,hy)
    ztaper = create_taper1d(tapz1,tapz2,nz,hz)
    ztaper = z * 0
    taper = np.zeros((nx,ny,nz))

    nh1 = int(np.floor(tapz1/hz))
    nh2 = int(np.floor(tapz2/hz))
    z1 = nh1+nh2-1
    ztaper[nh1:z1+1] = coswin(nh2,1)
    ztaper[z1+1:nz] = 1.

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                taper[i,j,k] = xtaper[i] * ytaper[j] * ztaper[k]
    
    intp = RegularGridInterpolator((x,y,z),taper,bounds_error=False,fill_value=0)

    # now loop every file to interpolate
    for myrank in range(nprocs):
        # read coordinates
        filename = './optimize/MODEL_M%02d'%(iter) + '/proc%06d'%(myrank) + '_external_mesh.bin'
        f = FortranFile(filename)
        nspec = f.read_ints('i4')[0]
        _ = f.read_ints('i4')[0]
        f.read_ints('i4')
        ibool = f.read_ints('i4') - 1
        xstore = f.read_reals('f4')
        ystore = f.read_reals('f4')
        zstore = f.read_reals('f4')
        f.close()

        # raed search direction 
        #out_list = ['dbulk','dbeta','drho']
        out_list = get_direc_name_list()
        for i in range(3):
            filename = './optimize/SUM_KERNELS_M%02d'%(iter)+ '/proc%06d'%(myrank)  \
                      + '_' + out_list[i] + '.bin'
            f = FortranFile(filename)
            kl = f.read_record('f4')
            kl *= intp((xstore[ibool],ystore[ibool],zstore[ibool]))

            kl  = np.float32(kl)
            filename = './optimize/SUM_KERNELS_M%02d'%(iter) + '/proc%06d'%(myrank)  \
                      + '_' + out_list[i] + '.bin'
            f = FortranFile(filename,"w")
            f.write_record(kl)
            f.close()

main()