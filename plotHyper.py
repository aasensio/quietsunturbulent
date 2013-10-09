import numpy as np
import matplotlib.pyplot as pl
import scipy.special as sp
import scipy.stats as stat
import scipy.integrate as integ
from matplotlib.ticker import MaxNLocator
import brewer2mpl

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


nTicks = 7

# Load the Markov chains that have been already thinned with thinning.py
f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

nLines = (npar-4) / 6

ch = np.load('parameters.npy')

quantiles = stat.mstats.mquantiles(ch[:,3*nLines:4*nLines], prob=[0.5-0.68/2, 0.5, 0.5+0.68/2], axis=0)

pl.close('all')

# Figure with the hyperparameters
fig2 = pl.figure(num=2, figsize=(14,9))
fig2.set_zorder(2)
pl.clf()

loop = 0
parameter = [r'$\gamma$',r'$\delta$']

nRows = 2
nCols = 3
a = np.arange(nRows*nCols)+1
order = a
#order = np.hstack(a.reshape(nRows, nCols).T)

for j in range(2):
	ax = fig2.add_subplot(nRows,nCols,order[loop])
	ax.plot(ch[:,-4+j], color='#969696')
	ax.xaxis.set_major_locator(MaxNLocator(nTicks))
	ax.set_xlabel('Iteration')
	ax.set_ylabel(parameter[j])
	loop += 1
		
maxB = 1000
nPoints = 1000

## Magnetic field strength
BMin = 1
logBMin = np.log10(BMin)
BMax = 1.2e3
logBMax = np.log10(BMax)
logB = np.linspace(logBMin, logBMax,nPoints)
B = 10.0**logB
alpha = ch[:,-4]
beta = ch[:,-3]
pB = LogNormalAvgPrior(B, alpha, beta)

ax = fig2.add_subplot(nRows,nCols,order[loop])
ax.semilogx(B,pB, color='#507FED')
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'p(B)')
#ax.xaxis.set_major_locator(MaxNLocator(nTicks))
ax.set_xlim(BMin,BMax)
loop += 1

pBCumulative = integ.cumtrapz(pB, B, initial=0)

ax = fig2.add_subplot(nRows,nCols,order[loop])
ax.semilogx(B,pBCumulative / pBCumulative.max(), color='#507FED')
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'cdf(B)')
ax.set_xlim(BMin,BMax)
#ax.xaxis.set_major_locator(MaxNLocator(nTicks))
ax.axhline(y=0.9,color='k',ls='dashed')
ax.axhline(y=0.95,color='k',ls='dashed')
loop += 1

cmap = pl.get_cmap('copper')

ax = fig2.add_subplot(nRows,nCols,order[loop])
for i in range(nLines):
	ax.semilogy(i, quantiles[1,i], 'o', color=cmap(float(i)/(nLines-1)))
	ax.errorbar(i, quantiles[1,i], yerr=[[quantiles[1,i]-quantiles[0,i]],[quantiles[2,i]-quantiles[1,i]]], color=cmap(float(i)/(nLines-1)))	

ax.set_ylim(BMin,BMax)
ax.set_xlabel('Line')
ax.set_ylabel('B [G]')
loop += 1

ax = fig2.add_subplot(nRows,nCols,order[loop])
ax.semilogy(pB,B, color='#507FED')
ax.set_ylim(BMin,BMax)
ax.set_ylabel('B [G]')
ax.set_xlabel('p(B)')
ax.xaxis.set_major_locator(MaxNLocator(4))

fig2.tight_layout()
fig2.savefig("hyperparameters.pdf")