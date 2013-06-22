#! /usr/bin/env python
############################################################
# Yigit Dallilar 13.06.2013                                #
# DTU Space - 25 perc modulation plot script               #
############################################################

import pylab as plt
import string as str
import numpy as np

f = open("15-165mod.txt").readlines()

a,b,c,d,e,ff = [],[],[],[],[],[]

for x in f:
    a.append(float(str.split(x)[0]))
    b.append(float(str.split(x)[1]))
    c.append(float(str.split(x)[2]))
    d.append(float(str.split(x)[3]))
    e.append(float(str.split(x)[4]))
    ff.append(float(str.split(x)[5]))

t = np.linspace(0,1,num=256)
a = np.array(a)
b = np.array(b)
c = np.array(c)
d = np.array(d)
e = np.array(e)
ff = np.array(ff)


plt.subplot(7,1,1)
plt.plot(t,a)
plt.grid(True)
plt.ylabel("A mod")
plt.subplot(7,1,2)
plt.plot(t,b)
plt.grid(True)
plt.ylabel("B mod")
plt.subplot(7,1,3)
plt.plot(t,c)
plt.grid(True)
plt.ylabel("C mod")
plt.subplot(7,1,4)
plt.plot(t,d)
plt.grid(True)
plt.ylabel("D mod")
plt.subplot(7,1,5)
plt.plot(t,e)
plt.grid(True)
plt.ylabel("E mod")
plt.subplot(7,1,6)
plt.plot(t,ff)
plt.grid(True)
plt.ylabel("F mod")
plt.subplot(7,1,7)
plt.plot(t,a+b+c+d+e+ff)
plt.grid(True)
plt.ylabel("SUM")
plt.xlabel("Theta (2*pi)")
plt.ylim([0,5])
plt.show()
