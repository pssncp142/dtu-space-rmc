#! /usr/bin/env python
############################################################
# Yigit Dallilar 13.06.2013                                #
# DTU Space - 25 perc modulation plot script               #
############################################################

import pylab as plt
import string as str
import numpy as np

f = open("25mod.txt").readlines()

a,b,c,d = [],[],[],[]

for x in f:
    a.append(float(str.split(x)[0]))
    b.append(float(str.split(x)[1]))
    c.append(float(str.split(x)[2]))
    d.append(float(str.split(x)[3]))

t = np.linspace(0,1,num=256)
a = np.array(a)
b = np.array(b)
c = np.array(c)
d = np.array(d)

plt.subplot(5,1,1)
plt.plot(t,a)
plt.grid(True)
plt.ylabel("A mod")
plt.subplot(5,1,2)
plt.plot(t,b)
plt.grid(True)
plt.ylabel("B mod")
plt.subplot(5,1,3)
plt.plot(t,c)
plt.grid(True)
plt.ylabel("C mod")
plt.subplot(5,1,4)
plt.plot(t,d)
plt.grid(True)
plt.ylabel("D mod")
plt.subplot(5,1,5)
plt.plot(t,a+b+c+d)
plt.grid(True)
plt.ylabel("SUM")
plt.xlabel("Theta (2*pi)")
plt.ylim([0,2])
plt.show()
