import numpy as np
import matplotlib.pyplot as pl
from scipy.special import wofz

"""
Generate an artificial atlas of Milne-Eddington lines with a certain magnetic field to see if
we are able to detect it.
"""

def voigt(x, y):
   # The Voigt function is also the real part of 
   # w(z) = exp(-z^2) erfc(iz), the complex probability function,
   # which is also known as the Faddeeva function. Scipy has 
   # implemented this function under the name wofz()
   z = x + 1j*y
   I = wofz(z).real
   return I

data = np.loadtxt('../lines.dat',skiprows=1)

lambdaAxis = np.asarray([])
spectrum = np.asarray([])

for i in range(data.shape[0]):
	lambda0 = data[i,2]
	if (lambda0 != 0):
		Gt = data[i,7]
		lambdaLeft = data[i,10] - 0.01
		lambdaRight = data[i,11] + 0.01
		
		deltaVDoppler = np.random.uniform(0.5,2.0)
		deltaLDoppler = lambda0 * deltaVDoppler * 1.e5 / 3.e10
		
		beta0 = np.random.uniform(0.5,20)
		etal = np.random.uniform(0.5,20)
		a = np.random.uniform(0.0,0.5)
		B = 500.0
		
		l = np.linspace(lambdaLeft, lambdaRight, 40)
			
		deltaT = np.sqrt(deltaLDoppler**2 + 4.0*Gt*(4.6686e-13*lambda0**2*B)**2)
				
		Psi = deltaLDoppler / (np.sqrt(np.pi)*deltaT) * voigt((l - lambda0)/deltaT, a*deltaLDoppler / deltaT)	
		line = (1.0+beta0 / (1.0 + etal*Psi)) / (1.0+beta0) + 0.01*np.random.randn(40)
		
		lambdaAxis = np.append(lambdaAxis, l)
		spectrum = np.append(spectrum, line)
	
pl.plot(spectrum)

f=open("atlasArtificial.bin","wb")

size=np.asarray(lambdaAxis.shape)
size.tofile(f)
res = lambdaAxis.astype(float)
res.tofile(f)
res = spectrum.astype(float)
res.tofile(f)
f.close()
