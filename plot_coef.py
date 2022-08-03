#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np


fp=open("coefs","r")

e_indxs=[]
e_vals=[]

a0_indxs=[]
a0_vals=[]
a1_indxs=[]
a1_vals=[]

while True:

    line=fp.readline()
    if not line:
        break

    strs=line.split()
    e_indxs.append(int(strs[0].strip()))
    e_vals.append(float(strs[1].strip()))

    line=fp.readline()
    if not line:
        break

    strs=line.split()
    a0_indxs.append(int(strs[0].strip()))
    a0_vals.append(float(strs[1].strip()))
    
    line=fp.readline()
    if not line:
        break

    strs=line.split()
    a1_indxs.append(int(strs[0].strip()))
    a1_vals.append(float(strs[1].strip()))


fig=plt.figure()
ax = fig.add_subplot(1,1,1)
fig1=plt.figure()
ax1 = fig1.add_subplot(1,1,1)

ax.plot(e_indxs)
ax.plot(a0_indxs)
ax.plot(a1_indxs)

ax1.plot(e_vals)
ax1.plot(a0_vals)
ax1.plot(a1_vals)


plt.show()
