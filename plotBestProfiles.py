import numpy as np
import matplotlib.pyplot as pl
import scipy.special as sp
from matplotlib.ticker import MaxNLocator


nTicks = 5

# Load the Markov chains that have been already thinned with thinning.py
f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

ch = np.load('parameters.npy')

pl.close('all')

# Lines
f, ax = pl.subplots(5, sharex=True, sharey=True, figsize=(12,12))
lines = np.loadtxt('best.profiles')
nLambda = lines[:,0].size
for i in range(5):
	left = i * nLambda / 5
	right = (i+1) * nLambda / 5
	ax[i].plot(lines[left:right,2], '.', color='#969696')
	ax[i].plot(lines[left:right,1], color= '#507FED')
	ax[i].set_ylim(0.1,1.05)
	ax[i].set_xlim(0,right-left)
f.subplots_adjust(hspace=0)
pl.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
big_ax = f.add_subplot(111)
big_ax.set_axis_bgcolor('none')
big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
big_ax.set_ylabel('Normalized intensity')
big_ax.set_xlabel('Wavelength point')
f.savefig("bestProfiles.pdf")