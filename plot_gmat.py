#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


f=open("hold","r")

ind0=[]
ind1=[]
val=[]
data=[]

line=f.readline()
while line:
    strs=line.split()

    if len(strs)==2:
        data.append(float(strs[1].strip()))
        
    if len(strs)==5:
        ind0.append(int(strs[0].strip()))
        ind1.append(int(strs[1].strip()))
        val.append(float(strs[4].strip()))
        
    line=f.readline()
        
neqn=len(data)
ind0=np.asarray(ind0)
ind1=np.asarray(ind1)
val=np.asarray(val)

ngrid=np.amax(ind1)+1
print(len(data),len(val),ngrid)

a_array=np.zeros((neqn,ngrid))
for j in range(len(val)):
    a_array[ind0[j],ind1[j]]=val[j]


b=np.matmul(a_array.T,a_array)
aa_i=np.linalg.pinv(b)

a_ai=np.matmul(aa_i,b)

print(a_array.shape)
print(aa_i.shape)
plt.contour(a_ai)
# for j in range(20):
#     plt.plot(a_ai[40*j+3300,:])

plt.show()
