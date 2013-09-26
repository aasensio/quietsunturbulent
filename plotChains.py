import numpy as np
import matplotlib.pyplot as pl
import scipy.special as sp
from matplotlib.ticker import MaxNLocator

def IGAvgPrior(x, alpha, beta):
	pf = np.zeros(len(x))
	for i in range(len(x)):		
		logy = alpha * np.log(beta) - sp.gammaln(alpha) - (alpha+1.0) * np.log(x[i]) - beta / x[i]
		pf[i] = np.mean(np.exp(logy))
	return pf


nTicks = 5

# Load the Markov chains that have been already thinned with thinning.py
f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

ch = np.load('parameters.npy')
ch = ch[-2000:,:]

pl.close('all')

# Figure with the chains
fig1 = pl.figure(num=1, figsize=(20,10))
pl.clf()

nTicks = 4
loop = 1
nCols = 6
nPars = 6
whichLines = [0,65,130,196,261,326]
labels = [r'$\beta_0$',r'$\eta_l$',r'$\Delta \lambda_D$',r'$B$',r'$a$',r'$\sigma$']
for i in range(nPars):	
	for j in range(nCols):
		ax = fig1.add_subplot(nPars,nCols,loop)
		ax.plot(ch[:,i*(npar-2)/nPars+whichLines[j]])
		ax.xaxis.set_major_locator(MaxNLocator(nTicks))
		if (j==0):
			ax.set_ylabel(labels[i])
		loop += 1
		
fig1.tight_layout()		
fig1.savefig("chains.pdf")