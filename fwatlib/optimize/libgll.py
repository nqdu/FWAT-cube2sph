from numba import jit 
import numpy as np 

@jit(nopython=True)
def legendre(n,x):
    leg = 1.
    if n == 0:
        leg = 1
    elif n==1:
        leg = x
    else:
        leg_down1 = x; leg_down2 = 1.
        for i in range(2,n+1):
            leg = (2*i-1)*x*leg_down1/i - (i-1)*leg_down2/i
            leg_down2 = leg_down1
            leg_down1 = leg

    return leg

@jit(nopython=True)
def dlegendre(n,x):
    dleg = 0.
    if n == 0:
        dleg = 0.
    elif n == 1:
        dleg = 1.
    else:
        leg_down1 = x; leg_down2 = 1.
        dleg_down1 = 1.; dleg_down2 = 0.
        for i in range(2,n+1):
            leg = (2*i-1)*x*leg_down1/i - (i-1)*leg_down2/i
            dleg = dleg_down2 + (2*i-1)*leg_down1
            leg_down2 = leg_down1
            leg_down1 = leg
            dleg_down2 = dleg_down1
            dleg_down1 = dleg

    return dleg

@jit(nopython=True)
def gauss_legendre_lobatto(length:int):
    #// define some constants
    tolerance = 4.0 * 1.0e-6
    nnewton_iter = 100
    pi = np.pi
    
    # allocate space
    n = length - 1; #// order of polynomial
    x = np.array([0.0 for i in range(length)])
    w = x * 1.
    if n == 1:
        x[0] = -1.; x[1] = 1.
        w[0] = 1.; w[1] = 1.
    else:
        leg = 0.; dleg = 0.; delta = 0.
        # set end points
        x[0]   = -1.0; x[n] = 1.
        w[0]   =  2./(n*(n+1.)); w[n] =  2./(n*(n+1.))

        for i in range(1,(n+1)//2):
            #// initial guess from an approximate form given by SV Parter (1999)
            x[i] = -np.cos( (i+0.25)*pi/n  - 3/(8*n*pi*(i+0.25)))

            #// newton iteration
            for j in range(nnewton_iter):
                leg = legendre(n+1,x[i]) - legendre(n-1,x[i])
                dleg = dlegendre(n+1,x[i]) - dlegendre(n-1,x[i])
                delta = -leg/dleg
                x[i] += delta
                if (np.abs(delta) <= tolerance * np.abs(x[i]) ): break
            
            x[n-i] = - x[i]
            leg = legendre(n, x[i])
            w[i] = 2./(n*(n+1.)*leg*leg)
            w[n-i] = w[i]

        if n %2 == 0 :
            x[n//2] = 0.
            leg = legendre(n, 0.0)
            w[n//2]  = 2./(n*(n+1.)*leg*leg)
    

    return x,w

# x,w = gauss_legendre_lobatto(5)
# print(x)
# print(w)
# np.savetxt("utils/libgll.txt",w,fmt='%f')

def get_gll_weights(N=5):
    if N == 5:
        w = np.array([
            0.100000,
            0.544444,
            0.711111,
            0.544444,
            0.100000])
    else:
        _,w = gauss_legendre_lobatto(N)

    return w
