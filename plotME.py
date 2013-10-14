import numpy as np
import matplotlib.pyplot as pl
import scipy.special as sp
import scipy.stats as stat
import scipy.integrate as integ
from matplotlib.ticker import MaxNLocator
import brewer2mpl

nTicks = 7

# Load the Markov chains that have been already thinned with thinning.py
f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

nLines = (npar-3) / 5

ch = np.load('parameters.npy')
nStep = ch.shape[0]

pl.close('all')
fig2 = pl.figure(num=2)
ax = fig2.add_subplot(1,1,1)
ax.plot(ch[:,-3])

