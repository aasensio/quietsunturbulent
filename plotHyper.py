import numpy as np
import matplotlib.pyplot as pl
import scipy.special as sp
import scipy.stats as stat
from matplotlib.ticker import MaxNLocator

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

def LogNormalAvgPrior(x, mu, sigma):
	pf = np.zeros(len(x))
	for i in range(len(x)):		
		logy = -np.log(sigma) - np.log(x[i]) - (np.log(x[i]) - mu)**2 / (2.0*sigma**2)
		pf[i] = np.mean(np.exp(logy))
	return pf


nTicks = 5

# Load the Markov chains that have been already thinned with thinning.py
f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

nLines = (npar-2) / 6

ch = np.load('parameters.npy')

quantiles = stat.mstats.mquantiles(ch[:,3*nLines:4*nLines], prob=[0.5-0.68/2, 0.5, 0.5+0.68/2], axis=0)

pl.close('all')

# Figure with the hyperparameters
fig2 = pl.figure(num=2, figsize=(10,6))
pl.clf()

loop = 1
parameter = [r'$\gamma$',r'$\delta$']

for j in range(2):
	ax = fig2.add_subplot(2,3,loop)
	ax.plot(ch[:,-2+j], color='#969696')
	ax.xaxis.set_major_locator(MaxNLocator(nTicks))
	ax.set_xlabel('Iteration')
	ax.set_ylabel(parameter[j])
	loop += 1
		
maxB = 1000
nPoints = 1000

## Magnetic field strength
logB = np.linspace(-3,3,nPoints)
B = 10**logB
alpha = ch[:,-2]
beta = ch[:,-1]
pB = LogNormalAvgPrior(B, alpha, beta)

ax = fig2.add_subplot(2,3,loop)
ax.plot(B,pB, color='#507FED')
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'p(B)')
ax.xaxis.set_major_locator(MaxNLocator(nTicks))
ax.set_xlim(0,maxB)
loop += 1

pBCumulative = np.cumsum(pB)

ax = fig2.add_subplot(2,3,loop)
ax.plot(B,pBCumulative / pBCumulative.max(), color='#507FED')
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'cdf(B)')
ax.set_xlim(0,maxB)
ax.xaxis.set_major_locator(MaxNLocator(nTicks))
ax.axhline(y=0.9,color='k',ls='dashed')
ax.axhline(y=0.95,color='k',ls='dashed')
loop += 1

ax = fig2.add_subplot(2,3,loop)
for i in range(nLines):
	ax.semilogy(i, quantiles[1,i], 'o')
	ax.errorbar(i, quantiles[1,i], yerr=[[quantiles[0,i]],[quantiles[2,i]]])

ax.set_ylim(1e0,1e3)
loop += 1

ax = fig2.add_subplot(2,3,loop)
ax.semilogy(pB,B, color='#507FED')
ax.set_ylim(1e0,1e3)

fig2.tight_layout()
fig2.savefig("hyperparameters.pdf")