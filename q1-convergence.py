import matplotlib.pyplot as pl
from numpy import genfromtxt

pl.clf()
pl.plot(-genfromtxt('const-100', usecols=0), label=r'$\alpha_k = 100$')
pl.plot(-genfromtxt('nonsum-100', usecols=0), label=r'$\alpha_k = 100 / \sqrt{k}$')
pl.plot(-genfromtxt('square-100', usecols=0), label=r'$\alpha_k = 100 / k$')
pl.legend(loc='lower right')
pl.xlabel('iteration')
pl.ylabel('profit')
pl.xlim(0,1000)
pl.savefig('q1-convergence.pdf')
