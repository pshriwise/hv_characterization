
from matplotlib import pyplot as plt
import numpy as np


#Retrieve data
raw = np.fromfile("data.dat",dtype=float,count=-1,sep="\t").reshape((2000,4))

#Average and remove duplicate values
raw = (raw[0::2]+raw[1::2])/2

#Get unique X,Y,Z values
X = np.unique(raw[:,0])

Y = np.unique(raw[:,1])

Z = raw[:,3].reshape(Y.shape[0],X.shape[0]) #third column is the split ratio value

Xv,Yv = np.meshgrid(X,Y)

fig, ax = plt.subplots(1,1)

ax.set_xlabel("Valance")
ax.set_ylabel("Area Fraction")

plt.pcolormesh(Xv,Yv,Z,)
plt.colorbar(format = '%2.0E')
plt.show()







