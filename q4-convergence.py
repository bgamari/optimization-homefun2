import matplotlib.pyplot as pl
import numpy as np

pl.plot(np.genfromtxt('u1-svt'), label='u1')
pl.plot(np.genfromtxt('u2-svt'), label='u2')
pl.plot(np.genfromtxt('u3-svt'), label='u3')
#pl.plot(np.genfromtxt('u4-svt'), label='u4')
pl.ylabel('relative error')
pl.xlabel('iteration')
pl.legend()
pl.savefig('q4-convergence.pdf')
