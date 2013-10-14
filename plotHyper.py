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

def LogNormalMixtureAvgPrior(x, mu1, sigma1, mu2, sigma2, p):
	pf = np.zeros(len(x))
	for i in range(len(x)):
		logp1 = np.log(p) - np.log(sigma1) - np.log(x[i]) - (np.log(x[i]) - mu1)**2 / (2.0*sigma1**2)
		logp2 = np.log(1.0-p) - np.log(sigma2) - np.log(x[i]) - (np.log(x[i]) - mu2)**2 / (2.0*sigma2**2)
		logp = np.logaddexp(logp1,logp2)
		pf[i] = np.mean(np.exp(logp))
	return pf


nTicks = 7

# Load the Markov chains that have been already thinned with thinning.py
f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

nLines = (npar-7) / 6

ch = np.load('parameters.npy')
nStep = ch.shape[0]

p = ch[:,-1]
sigma2 = ch[:,-2]
mu2 = ch[:,-3]
sigma1 = ch[:,-6]
mu1 = ch[:,-7]

hyperp = [mu1, sigma1, mu2, sigma2, p]

quantiles = stat.mstats.mquantiles(ch[:,3*nLines:4*nLines], prob=[0.5-0.68/2, 0.5, 0.5+0.68/2], axis=0)

pl.close('all')

# Figure with the hyperparameters
fig2 = pl.figure(num=2, figsize=(14,9))
fig2.set_zorder(2)
pl.clf()

loop = 0
parameter = [r'$\mu_1$',r'$\sigma_1$',r'$\mu_2$',r'$\sigma_2$','p']

nRows = 2
nCols = 3
a = np.arange(nRows*nCols)+1
order = a
#order = np.hstack(a.reshape(nRows, nCols).T)

ax = fig2.add_subplot(nRows,nCols,order[loop])
ax.plot(mu1, color='#969696')
ax.plot(mu2, color='#507FED')
ax.set_xlim(0,nStep)
ax.set_xlabel('Iteration')
ax.set_ylabel(r'$\mu$')
loop += 1

ax = fig2.add_subplot(nRows,nCols,order[loop])
ax.plot(sigma1, color='#969696')
ax.plot(sigma2, color='#507FED')
ax.set_xlim(0,nStep)
ax.set_xlabel('Iteration')
ax.set_ylabel(r'$\gamma$')
loop += 1

ax = fig2.add_subplot(nRows,nCols,order[loop])
ax.plot(p, color='#969696')
ax.plot(1.0-p, color='#507FED')
ax.set_xlim(0,nStep)
ax.set_xlabel('Iteration')
ax.set_ylabel('p')
loop += 1



nPoints = 1000

# Magnetic field strength
BMin = 1e-4
BMinPlot = 1.e-4
logBMin = np.log10(BMin)
BMax = 1.2e4
BMaxPlot = 1.2e4
logBMax = np.log10(BMax)
logB = np.linspace(logBMin, logBMax,nPoints)
B = 10.0**logB
pB = LogNormalMixtureAvgPrior(B, mu1, sigma1, mu2, sigma2, p)

ax = fig2.add_subplot(nRows,nCols,order[loop])
ax.loglog(B,pB, color='#507FED')
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'p(B)')
ax.set_xlim(BMinPlot,BMaxPlot)
loop += 1

pBCumulative = integ.cumtrapz(pB, B, initial=0)

ax = fig2.add_subplot(nRows,nCols,order[loop])
ax.set_xlim(BMinPlot,BMaxPlot)
ax.semilogx(B,pBCumulative / pBCumulative.max(), color='#507FED')
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'cdf(B)')
#ax.xaxis.set_major_locator(MaxNLocator(nTicks))
values = [np.interp(quantile,pBCumulative / pBCumulative.max(),B) for quantile in [0.5,0.9,0.1]]
ax.axhline(y=0.5,color='k',ls='dashed')
ax.text(2e-4, 0.52, '50% - '+"{:.0f}".format(values[0])+' G')
ax.axhline(y=0.9,color='k',ls='dashed')
ax.text(2e-4, 0.92, '90% - '+"{:.0f}".format(values[1])+' G')
ax.axhline(y=0.1,color='k',ls='dashed')
ax.text(2e-4, 0.12, '10% - '+"{:.0f}".format(values[2])+' G')
loop += 1

cmap = pl.get_cmap('copper')

ax = fig2.add_subplot(nRows,nCols,order[loop])
for i in range(nLines):
	ax.semilogy(i, quantiles[1,i], 'o', color=cmap(float(i)/(nLines-1)))
	ax.errorbar(i, quantiles[1,i], yerr=[[quantiles[1,i]-quantiles[0,i]],[quantiles[2,i]-quantiles[1,i]]], color=cmap(float(i)/(nLines-1)))	

ax.set_ylim(BMinPlot,BMaxPlot)
ax.set_xlabel('Line')
ax.set_ylabel('B [G]')

fig2.tight_layout()
#fig2.savefig("hyperparameters.pdf")