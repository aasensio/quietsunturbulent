import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from matplotlib.ticker import MaxNLocator

def IGAvgPrior(x, alpha, beta):
	pf = np.zeros(len(x))
	for i in range(len(x)):		
		logy = alpha * np.log(beta) - sp.gammaln(alpha) - (alpha+1.0) * np.log(x[i]) - beta / x[i]
		pf[i] = np.mean(np.exp(logy))
	return pf


nTicks = 5

f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

ch = np.fromfile('test.extract', dtype=np.float64).reshape((npar,nstep), order="FORTRAN")

length = len(ch[0,:])
whichElements = np.arange(0, length, 5)
ch = ch[:,whichElements]
ch = ch[:,-2000:]

plt.close('all')

# Figure with the chains
fig1 = plt.figure(num=1, figsize=(20,10))
plt.clf()

nTicks = 4
loop = 1
nCols = 7
nPars = 6
labels = [r'$\beta_0$',r'$\eta_l$',r'$\Delta \lambda_D$',r'$B$',r'$a$',r'$\sigma$']
for i in range(nPars):	
	for j in range(nCols):
		ax = fig1.add_subplot(nPars,nCols,loop)
		ax.plot(ch[i*(npar-2)/nPars+j,:])
		ax.xaxis.set_major_locator(MaxNLocator(nTicks))
		ax.set_ylabel(labels[i])
		loop += 1
		
fig1.tight_layout()		

# Figure with the hyperparameters
fig2 = plt.figure(num=2)
plt.clf()

loop = 1
for j in range(2):
	ax = fig2.add_subplot(2,2,loop)
	ax.plot(ch[-2+j,:])
	ax.xaxis.set_major_locator(MaxNLocator(nTicks))
	loop += 1
		
## Magnetic field strength
logB = np.linspace(-3,4,1000)
B = 10**logB
pB = np.zeros(100)
alpha = ch[-2,:]
beta = ch[-1,:]
pB = IGAvgPrior(B, alpha, beta)

ax = fig2.add_subplot(2,2,loop)
ax.plot(B,pB)
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'p(B)')
ax.xaxis.set_major_locator(MaxNLocator(nTicks))
ax.set_xlim(0,1000)
loop += 1

pBCumulative = np.cumsum(pB)

ax = fig2.add_subplot(2,2,loop)
ax.plot(B,pBCumulative / pBCumulative.max())
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'cdf(B)')
ax.set_xlim(0,250)
ax.xaxis.set_major_locator(MaxNLocator(nTicks))
ax.axhline(y=0.9,color='k',ls='dashed')
ax.axhline(y=0.95,color='k',ls='dashed')

# Lines
fig3 = plt.figure(num=3)
lines = np.loadtxt('best.profiles')
ax = fig3.add_subplot(1,1,1)
ax.plot(lines[:,2], 'r.')
ax.plot(lines[:,1])

