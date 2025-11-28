#!/bin/bash 

set -e
source parameters.sh 
source module_env_gmt

bounds=-R$LON0_H/$LON1_H/$LAT0_H/$LAT1_H
proj=-JM12C

if [ $# -ne 1 ]; then
    echo "Usage: get_mask.sh do_mask (=0/1)"
    exit 1
fi

do_mask=$1

python << EOF
import numpy as np 
n = 256
lat = np.linspace($LAT0_H,$LAT1_H,n)
lon = np.linspace($LON0_H,$LON1_H,n)

f = open("points.xyz", "w")
for i in range(n):
    for j in range(n):
        f.write("%f %f %d\n" % (lon[i], lat[j],i*n+j))
f.close()
EOF

# selection
if [ $do_mask -ne 0 ]; then 
  echo "Select points inside the region ..."
  gmt select points.xyz $bounds -Dh -Nk/s  > out.txt 
else 
  cat points.xyz > out.txt
fi

# read back
python << EOF
import numpy as np
data = np.loadtxt("out.txt")
idx = np.asarray(data[:,2], dtype=int)
points = np.loadtxt("points.xyz")
points[:,-1] = 1.
points[idx,-1] = np.nan

np.savetxt("profiles/mask.dat", points)
EOF

rm out.txt points.xyz