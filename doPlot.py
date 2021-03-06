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
#ch = ch[-2000:,:]

pl.close('all')

# Figure with the chains
fig1 = pl.figure(num=1, figsize=(20,10))
pl.clf()

nTicks = 4
loop = 1
nCols = 7
nPars = 6
labels = [r'$\beta_0$',r'$\eta_l$',r'$\Delta \lambda_D$',r'$B$',r'$a$',r'$\sigma$']
for i in range(nPars):	
	for j in range(nCols):
		ax = fig1.add_subplot(nPars,nCols,loop)
		ax.plot(ch[:,i*(npar-2)/nPars+j+300])
		ax.xaxis.set_major_locator(MaxNLocator(nTicks))
		ax.set_ylabel(labels[i])
		loop += 1
		
fig1.tight_layout()		

# Figure with the hyperparameters
fig2 = pl.figure(num=2)
pl.clf()

loop = 1
for j in range(2):
	ax = fig2.add_subplot(2,2,loop)
	ax.plot(ch[:,-2+j])
	ax.xaxis.set_major_locator(MaxNLocator(nTicks))
	loop += 1
		
## Magnetic field strength
logB = np.linspace(-3,1,1000)
B = 10**logB
pB = np.zeros(100)
alpha = ch[:,-2]
beta = ch[:,-1]
pB = IGAvgPrior(B, alpha, beta)

ax = fig2.add_subplot(2,2,loop)
ax.plot(B,pB)
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'p(B)')
ax.xaxis.set_major_locator(MaxNLocator(nTicks))
ax.set_xlim(0,0.2)
loop += 1

pBCumulative = np.cumsum(pB)

ax = fig2.add_subplot(2,2,loop)
ax.plot(B,pBCumulative / pBCumulative.max())
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'cdf(B)')
ax.set_xlim(0,0.2)
ax.xaxis.set_major_locator(MaxNLocator(nTicks))
ax.axhline(y=0.9,color='k',ls='dashed')
ax.axhline(y=0.95,color='k',ls='dashed')

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