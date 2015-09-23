
import sys
from matplotlib import pyplot as plt
import numpy as np
import argparse as ap

parser = ap.ArgumentParser()

parser.add_argument('--datafile', type=str, default='data.dat', help='Datafile containing timing information.')
parser.add_argument('--title',dest='plt_title', type=str,
                    default='High Valance Characterizaion Timing Plot',
                    help = 'Used to set the title of the created plot')

args = parser.parse_args(sys.argv[1:])


#Retrieve data
raw = np.fromfile(args.datafile,dtype=float,count=-1,sep="\t")
raw = raw.reshape((raw.shape[0]/4,4))

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

plt.title(args.plt_title)
plt.pcolormesh(Xv,Yv,Z,)
plt.colorbar(format = '%2.0E')
plt.show()







