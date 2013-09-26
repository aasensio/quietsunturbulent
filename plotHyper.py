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

# Figure with the hyperparameters
fig2 = pl.figure(num=2)
pl.clf()

loop = 1
parameter = [r'$\gamma$',r'$\delta$']

for j in range(2):
	ax = fig2.add_subplot(2,2,loop)
	ax.plot(ch[:,-2+j], color='#969696')
	ax.xaxis.set_major_locator(MaxNLocator(nTicks))
	ax.set_xlabel('Iteration')
	ax.set_ylabel(parameter[j])
	loop += 1
		
## Magnetic field strength
logB = np.linspace(-3,1,1000)
B = 10**logB
pB = np.zeros(100)
alpha = ch[:,-2]
beta = ch[:,-1]
pB = IGAvgPrior(B, alpha, beta)

ax = fig2.add_subplot(2,2,loop)
ax.plot(B,pB, color='#507FED')
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'p(B)')
ax.xaxis.set_major_locator(MaxNLocator(nTicks))
ax.set_xlim(0,0.2)
loop += 1

pBCumulative = np.cumsum(pB)

ax = fig2.add_subplot(2,2,loop)
ax.plot(B,pBCumulative / pBCumulative.max(), color='#507FED')
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'cdf(B)')
ax.set_xlim(0,0.08)
ax.xaxis.set_major_locator(MaxNLocator(nTicks))
ax.axhline(y=0.9,color='k',ls='dashed')
ax.axhline(y=0.95,color='k',ls='dashed')
fig2.tight_layout()
fig2.savefig("hyperparameters.pdf")