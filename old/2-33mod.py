#! /usr/bin/env python
############################################################
# Yigit Dallilar 13.06.2013                                #
# DTU Space - 33 perc modulation plot script               #
############################################################

import pylab as plt
import string as str
import numpy as np

f = open("2-33mod.txt").readlines()

a,b,c = [],[],[]

for x in f:
    a.append(float(str.split(x)[0]))
    b.append(float(str.split(x)[1]))
    c.append(float(str.split(x)[2]))

t = np.linspace(0,1,num=256)
a = np.array(a)
b = np.array(b)
c = np.array(c)

plt.subplot(4,1,1)
plt.plot(t,a)
plt.grid(True)
plt.ylabel("A mod")
plt.subplot(4,1,2)
plt.plot(t,b)
plt.grid(True)
plt.ylabel("B mod")
plt.subplot(4,1,3)
plt.plot(t,c)
plt.grid(True)
plt.ylabel("C mod")
plt.subplot(4,1,4)
plt.plot(t,a+b+c)
plt.grid(True)
plt.ylabel("SUM")
plt.xlabel("Theta (2*pi)")
plt.ylim([0,5])
plt.show()
