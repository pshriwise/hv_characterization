
import sys
from matplotlib import pyplot as plt
from matplotlib import colors
import numpy as np
import argparse as ap

parser = ap.ArgumentParser()

parser.add_argument('--datafile', type=str, default='data.dat', help='Datafile containing timing information.')
parser.add_argument('--title',dest='plt_title', type=str,
                    default='High Valance Characterizaion Timing Plot',
                    help = 'Used to set the title of the created plot')
parser.add_argument('--no-log',dest='no_log',action='store_true',default=False,help='Indicates that the colorbar should not be in log scale.')

args = parser.parse_args(sys.argv[1:])

#Retrieve data from specified datafile
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

if args.no_log:
    mesh = plt.pcolormesh(Xv,Yv,Z)
else:
    mesh = plt.pcolormesh(Xv,Yv,Z,norm=colors.LogNorm())

mesh.set_clim(1e-6,1e-03)
cb = plt.colorbar(format = '%2.0E')
cb.set_label("Time (s)")
plt.show()







