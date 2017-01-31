
import pandas as pd
import numpy as np
import plotly.offline as py
import plotly.graph_objs as go
import plotly
import math
plotly.offline.init_notebook_mode()

d = pd.read_csv('./bin/data.dat',delim_whitespace=True)

#Retrieve data from specified datafile
raw = np.fromfile('./bin/data.dat',dtype=float,count=-1,sep="\t")
raw = raw.reshape((raw.shape[0]/4,4))

#Average and remove duplicate values
#raw = (raw[0::2]+raw[1::2])/2

#Get unique X,Y,Z values
X = np.unique(raw[:,0])

Y = np.unique(raw[:,1])

Z = raw[::2,3].reshape(Y.shape[0],X.shape[0]) #third column is the split ratio value

print Z

data = [go.Heatmap( z=Z, colorscale='Viridis')]


data = [{
    'z': 1.e6*Z,
    'type': 'heatmap',
    'colorscale': [
        [0, 'rgb(250, 250, 250)'],        #0
        [1./10000, 'rgb(200, 200, 200)'], #10
        [1./1000, 'rgb(150, 150, 150)'],  #100
        [1./100, 'rgb(100, 100, 100)'],   #1000
        [1./10, 'rgb(50, 50, 50)'],       #10000
        [1., 'rgb(0, 0, 0)'],             #100000
    ],
    'colorbar': {
        'tick0': 0,
        'tickmode': 'array',
        'tickvals': [0.001, 0.0001, 0.00001, 0.000001]
    }
    }]

fig = {'data': data}

py.plot(data)
