#! /usr/bin/env python
############################################################
# Yigit Dallilar 13.06.2013                                #
# DTU Space - 33 perc modulation plot script               #
############################################################

import pylab as plt
import string as str
import numpy as np

f = open("50mod.txt").readlines()

a,b,c = [],[],[]

for x in f:
    a.append(float(str.split(x)[0]))
    b.append(float(str.split(x)[1]))

t = np.linspace(0,1,num=256)
a = np.array(a)
b = np.array(b)

plt.subplot(3,1,1)
plt.plot(t,a)
plt.grid(True)
plt.ylabel("A mod")
plt.subplot(3,1,2)
plt.plot(t,b)
plt.grid(True)
plt.ylabel("B mod")
plt.subplot(3,1,3)
plt.plot(t,a+b)
plt.grid(True)
plt.ylabel("SUM")
plt.xlabel("Theta (2*pi)")
plt.ylim([0,2])
plt.show()
