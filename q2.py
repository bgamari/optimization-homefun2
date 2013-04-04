#!/usr/bin/python

import matplotlib.pyplot as pl
import numpy as np

h = np.genfromtxt('h', usecols=(0,1))
h2 = np.genfromtxt('h2')

n = 80
h = h[:n,:]
h2 = h2[:n,:]
g = np.empty((2*n,2))
g[0::2,:] = h
g[1::2,:] = h2

obj = lambda x1,x2: x1 + 2*x2 - 1

pl.subplot(211)
pl.plot(g[:,0], g[:,1])
pl.plot([0,1], [1,0], c='k')
pl.xlabel('$x_1$')
pl.ylabel('$x_2$')

pl.subplot(212)
pl.plot(obj(g[:,0], g[:,1]))
pl.xlabel('iteration')
pl.ylabel(r'$f(\vec x)$')
pl.savefig('q2-convergence.pdf')
