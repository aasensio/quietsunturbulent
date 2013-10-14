import numpy as np
import matplotlib.pyplot as pl
import scipy.special as sp
from matplotlib.ticker import MaxNLocator
import scipy.stats as stat

def IGAvgPrior(x, alpha, beta):
	pf = np.zeros(len(x))
	for i in range(len(x)):		
		logy = alpha * np.log(beta) - sp.gammaln(alpha) - (alpha+1.0) * np.log(x[i]) - beta / x[i]
		pf[i] = np.mean(np.exp(logy))
	return pf

def GaussAvgPrior(x, alpha, beta):
	pf = np.zeros(len(x))
	for i in range(len(x)):		
		logy = -(x[i]-alpha)**2 / (2.0*beta**2) - np.log(beta)
		pf[i] = np.mean(np.exp(logy))
	return pf


nTicks = 5

# Load the Markov chains that have been already thinned with thinning.py
f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())
nLines = (npar-7) / 6

ch = np.load('parameters.npy')

pl.close('all')

# Figure with the chains
fig1 = pl.figure(num=1, figsize=(20,10))
pl.clf()

nTicks = 4
loop = 1
nCols = 7
nPars = 6
whichLines = [0,1,2,3,4,5]
xPos = np.arange(nLines)
labels = [r'$\beta_0$',r'$\eta_l$',r'$\Delta v_D$',r'$B$',r'$a$',r'$\sigma$']
for i in range(nPars):	
	for j in range(nCols-1):
		ax = fig1.add_subplot(nPars,nCols,loop)
		ax.plot(ch[:,i*(npar-2)/nPars+whichLines[j]])
		ax.xaxis.set_major_locator(MaxNLocator(nTicks))
		if (j==0):
			ax.set_ylabel(labels[i])
		loop += 1
	ax = fig1.add_subplot(nPars,nCols,loop)
	quantiles = stat.mstats.mquantiles(ch[:,i*nLines:(i+1)*nLines], prob=[0.5-0.68/2, 0.5, 0.5+0.68/2], axis=0)
	for j in range(nLines):
		ax.plot(j, quantiles[1,j], 'o')
		ax.errorbar(j, quantiles[1,j], yerr=[[quantiles[1,i]-quantiles[0,i]],[quantiles[2,i]-quantiles[1,i]]])
	loop += 1
		
fig1.tight_layout()		
fig1.savefig("chains.pdf")