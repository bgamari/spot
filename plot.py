#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as pl

a = np.genfromtxt('h', dtype=None, names='r,phi,z,i')
x = a['r'] * np.cos(a['phi'])
y = a['r'] * np.sin(a['phi'])
pl.tripcolor(x, y, a['i'])
pl.show()
