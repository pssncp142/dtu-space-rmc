#! /usr/bin/env python                                         
################################################################
# Yigit Dallilar 08.06.2013                                    #
# Plotting script                                              #
################################################################

import pylab as plt
import string as str
import numpy as np

a,b=[],[]
tmp = open("realdata.txt").readlines()

for x in tmp :
    a.append(float(str.split(x)[0]))
    b.append(float(str.split(x)[1]))

a = np.array(a)
b = np.array(b)
c = np.linspace(0.,1.,256)

fig = plt.figure()
fig.add_subplot(311)
plt.plot(c,a)
plt.grid(True)
plt.ylabel("Count A")
plt.ylim([0.,60.])
fig.add_subplot(312)
plt.plot(c,b)
plt.grid(True)
plt.ylabel("Count B")
plt.ylim([0.,60.])
fig.add_subplot(313)
plt.plot(c,a+b)
plt.ylabel("Count A & B")
plt.grid(True)
plt.ylim([0.,60.])
plt.xlabel("Theta (2*pi)")
plt.show()

################################################################
