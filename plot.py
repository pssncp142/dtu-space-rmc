#! /usr/bin/env python                                         
################################################################
# Yigit Dallilar 08.06.2013                                    #
# Plotting script                                              #
################################################################

import pylab as plt
import string as str
import numpy as np

a=[]
tot = np.zeros(256)
tmp = open("strip.txt").readlines()

sp = len(str.split(tmp[0]))

for x in tmp :
    for y in str.split(x[:-1]) :
        a.append(float(y))

a = np.array(a)
b = np.linspace(0.,1.,num=256)

i=0
while(i<sp) :
    plt.subplot(sp+1,1,i+1)
    plt.plot(b,a[i:sp*256:sp])
    plt.grid(True)
    tot = tot + a[i:sp*256:sp]
    i = i + 1

plt.subplot(sp+1,1,sp+1)
plt.plot(b,tot)
#plt.ylim(0,4)
plt.grid(True)
plt.show()

################################################################
