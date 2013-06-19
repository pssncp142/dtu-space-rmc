#! /usr/bin/env python                                         
################################################################
# Yigit Dallilar 08.06.2013                                    #
# Plotting script                                              #
################################################################

import pylab as plt
import string as str
import numpy as np

a,b,c=[],[],[]
tmp = open("33real.txt").readlines()

for x in tmp :
    a.append(float(str.split(x)[0]))
    b.append(float(str.split(x)[1]))
    c.append(float(str.split(x)[2]))

a = np.array(a)
b = np.array(b)
c = np.array(c)
t = np.linspace(0.,1.,num=256)

plt.subplot(4,1,1)
plt.plot(t,a)
plt.grid(True)
plt.ylabel("Count A")
plt.subplot(4,1,2)
plt.plot(t,b)
plt.grid(True)
plt.ylabel("Count B")
plt.subplot(4,1,3)
plt.plot(t,c)
plt.grid(True)
plt.ylabel("Count C")
plt.subplot(4,1,4)
plt.plot(t,a+b+c)
plt.ylabel("SUM")
plt.grid(True)
plt.xlabel("Theta (2*pi)")
plt.show()

################################################################
