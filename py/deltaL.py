import numpy as np
import matplotlib.pyplot as pl

l = np.linspace(4300,7000,200)

for i in range(200):
	T = np.random.uniform(3000,7000,100)
	vmic = np.random.uniform(0,3e5,100)
	v = np.sqrt(2.0 * 1.381e-16 * T / (56.0 * 1.66e-24))
	x = np.tile(l[i],100)
	pl.plot(x, v*l[i]/3e10, '.')