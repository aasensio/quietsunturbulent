import numpy as np

f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

memMap = np.memmap('test.extract',dtype=np.float64,mode='r',shape=(nstep,npar))
ch = memMap[::4,:]

np.save('parameters.npy', ch)