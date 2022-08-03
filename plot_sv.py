#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


f=open("singVals","r")

npts=int(f.readline().strip())

svals=f.readline()
svals=svals.split()

svs=[]
for sv in svals:
    print(sv)
    svs.append(float(sv.strip()))

svs=np.asarray(svs)

plt.plot(svs)
plt.yscale("log")
plt.show()
