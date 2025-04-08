import numpy as np 
from libpca import PCA as PCA1 
import matplotlib.pyplot as plt 
from sklearn.decomposition import PCA

def do_pca(rec):
    # rec1 shape (nrec,npts)
    rec1 = rec - np.mean(rec,axis=1,keepdims=True)
    nt = rec1.shape[1]

    # eg,_ = np.linalg.eig(rec1.T @ rec)
    # print(eg[:5])

    # pca = PCA(n_components=rec.shape[0])
    # pca.fit(rec1)
    # stf = pca.components_[0,:]
    # print(pca.explained_variance_[:5]**2)
    # # u,s,v = np.linalg.svd(rec1)

    w,v = np.linalg.eig(rec1 @ rec1.T / nt)
    print(w[:5].real)
    
    # rec_new = u[:,1]
    stf = rec1.T @ v[:,0]
    
    # get flag
    sig1 = np.sign(rec1.flatten()[np.argmax(abs(rec1.flatten()))])
    sig2 = np.sign(stf[np.argmax(abs(stf))])
    stf *= sig1 * sig2

    return stf.real

def do_pca_avg(rec):
    stf1 = np.mean(rec,axis=0,keepdims=True)
    delta = np.max(abs(rec)) * 0.01

    for i in range(50):
        a = np.sum((stf1 - rec))**2
        grad = a / (np.sqrt(1 + a**2 / delta**2)) * np.sum(stf1 - rec,axis=0,keepdims=True)
        stf1 -= 0.01 * grad

    return stf1[0,:]

n = 100
nt = 1000
rec = np.zeros((n,nt))
t= np.linspace(0,1,nt)
for i in range(n):
    rec[i,:] = np.exp(-(t-0.5)**2/0.005) +   \
                0.1 * np.exp(-(t-0.5 + 0.2 *i / n)**2 / 0.005) +   \
                0.1 * np.random.randn(nt)

rec = rec - np.mean(rec,axis=1,keepdims=True)
stf1 = PCA1(np.float32(rec))
stf = do_pca_avg(rec)
print(np.dot(stf,stf),np.dot(stf1,stf1))


plt.figure(1)
plt.plot(t,stf1,label='f')
plt.plot(t,stf + 1,label='py')
plt.plot(t,np.mean(rec,axis=0),label='avg')
# plt.plot(t,rec[0,:])
plt.legend()
plt.savefig("test.jpg")
