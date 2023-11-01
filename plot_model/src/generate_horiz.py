import numpy as np
import sys  

def main():
    if len(sys.argv)!=9:
        print("Usage ./this lonmin lonmax latmin latmax z nlon nlat outfile")
        exit(1)

    lonmin,lonmax,latmin,latmax,z = map(lambda x:float(x),sys.argv[1:6])
    nlon,nlat = map(lambda x:int(x),sys.argv[6:8])
    outfile = sys.argv[8]
    lon = np.linspace(lonmin,lonmax,nlon)
    lat = np.linspace(latmin,latmax,nlat)

    f = open(outfile,"w")
    for i in range(nlon):
        for j in range(nlat):
            f.write("%f %f %f\n"%(lon[i],lat[j],z*1000))

    f.close()

main()

