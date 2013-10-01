import numpy as np

f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

memMap = np.memmap('test.extract',dtype=np.float64,mode='r',shape=(nstep,npar))
ch = memMap[-nstep/2::4,:]

np.save('parameters.npy', ch)