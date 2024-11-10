import numpy as np 
def determin_ntsamp(KERNEL_SPP,KERNEL_T0,dt):
    #! determine NSTEP_PER_FORWARD_OUTPUT based on the two parameters
    sit = np.floor(KERNEL_T0 / KERNEL_SPP / dt)
    dt_c = max(dt,dt*sit)
    samp = int(KERNEL_T0 / dt_c)
    sit = np.floor(KERNEL_T0 / samp / dt)
    
    return sit
  
print(determin_ntsamp(8,10.0,0.05))